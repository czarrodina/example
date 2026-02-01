import os
import argparse
import csv
import sys
import json
from ete3 import NCBITaxa
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from Properties import Properties, DEF_PROPERTIES_FILE

def main():
    # this program dedicated to identifying region(s) of a viral assembly that are dissimilar to a given taxa of interest
    # it does this by determining if blast hits for a segment are within a given rank and homology score to the taxa of interest
    # the program determines whether or not keep taxa within a given segment by determining how close a taxa is to the highest hit within the segment (-p percent)
    # the fuzzy top hit allows for both accounting for the low number of viruses sequenced relative to the total number of viruses,
    # the relatively high mutation rates/recombination rates of some viruses, and the relatively high amount of ambiguous bases in some viral assemblies
    # a high quality threshold is also used (-s score), which is a cutoff the represents a discrete threshold for scoring and for keeping segments
    # segments with high quality hits, that also pass the percent threshold, to taxa within the same rank cutoff are not saved
    # regions that are dissimilar are saved, combined with overlapping regions, and then scored to give a confidence value (with coverage/quality issues impacting score)
    # deletions relative to a reference and segments without blast hits (that are not low complexity) are also scored
    # the confidence value is an estimate as to likely a segment is to not be natural in origin (ie. synthetic), with 1 being lowest saved and 10 being highest
    parser = argparse.ArgumentParser(
        description="Identifies regions of a virus that have low homology to known isolates of the virus")
    parser.add_argument('-t', '--taxaid', dest='taxaid', type=int, required=True, help="Taxa id for the virus")
    parser.add_argument('-s', '--score', dest='threshold', type=float, default=0.8,
                        help="Score threshold for high quality hit")
    parser.add_argument('-i', '--infilename', dest='infilename', required=True, help="Name/Location of the input file")
    parser.add_argument('-o', '--outfilename', dest='outfilename', required=True,
                        help="Name/Location of the output file")
    parser.add_argument('-r', '--rank', dest='rank', type=str, default="genus", help="Rank cutoff in treating a taxa as the same as the taxa of interest")
    parser.add_argument('-p', '--percent', dest='percent', type=float, default=0.06,
                        help='Cutoff for percent difference between the highest score and the lowest kept score')
    parser.add_argument('-f', '--fasta', dest='fasta', required=True, help="Name/Location of the fasta file")
    parser.add_argument('-c', '--coverage', dest='coverage', required=False, default=None,
                        help="Name/Location of the low coverage file")
    parser.add_argument('-u', '--unmapped', dest='unmapped', required=False, default=None,
                        help="Name/Location of the unmapped segments file")
    parser.add_argument('-d', '--deletion', dest='deletion', required=False, default=None,
                        help="Name/Location of the deleted segments file")
    parser.add_argument('-n', '--figurename', dest='figname', required=False, default=None,
                        help="Leading name for the Figures")
    parser.add_argument('-e', '--etedata', dest='etedata', required=False, default=None,
                        help="Location of the ete NCBI database")
    parser.add_argument('--docker', dest="docker",
                        help="execution is going to take place inside a docker/singularity container.",
                        action='store_true')

    args = parser.parse_args()
    args_map = vars(args)
    def_properties_path = os.path.join(Properties.src_dir(), DEF_PROPERTIES_FILE)
    properties = Properties(def_properties_path, args_map)

    # set up NCBI taxonomy and find lineage for taxa of interest
    ete_db_path = args.etedata if args.etedata else properties.taxdb()
    ncbi = NCBITaxa(ete_db_path)
    lineage = ncbi.get_lineage(args.taxaid)
    errfile = open("log.err", 'w')
    # taxa id cutoff to be treated as the same taxa
    cutoff = 1
    # taxa cutoff rank of interest found
    run = False

    # seachers lineage for taxa rank cutoff
    for key in lineage:
        nrank = ncbi.get_rank([int(key)])
        if nrank[key] == args.rank:
            cutoff = key
            run = True
            ranklin = ncbi.get_lineage(cutoff)

    # if taxa rank cutoff found, run the program; otherwise print error that rank was not found
    if run == True:
        # saves cutoff properties for the given run
        properties = []
        properties.append(args.rank)
        properties.append(args.threshold)
        properties.append(args.percent)
        # array containing all taxa and scores for a given segment
        segment = []
        # array containing all taxa and scores for segments that are not within rank of taxa of interest
        data = []
        # array of segment scores 1-10 representing confidence that any given saved segment is worth looking at for synthetic alterations
        scores = []
        # justification for why the segment scores were chosen
        justification = []
        # data used to generate any figures (if figure file is given)
        graphdata = {}
        # number of kept segments
        count = 0
        # maximum homology score in a given segment
        maximum = 0.0
        # array containing anomolies in coverage (if file is given)
        coveragedata = None
        # json formatted data containing all segments of interest
        virus_data = {}
        # sequence for the assembly
        fastafile = open(args.fasta, 'r')
        sequencelist = parsefasta(fastafile)

        # parses coverage file to save any anomolies
        if args.coverage:
            coveragefile = open(args.coverage, 'r')
            coveragedata = parsecoverage(coveragefile)
            coveragefile.close()

        # file containing blast output for each segment
        csvfile = open(args.infilename, 'r')
        # json formatted file containing all segments of interest
        outfile = open(args.outfilename, 'w')

        # reads csv file and saves location important fields
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
        pident = fields.index("pident")
        squal = fields.index("subj_qual")
        status = fields.index("status")
        source = fields.index("source_taxid")
        dist = fields.index("distance_to_common")
        title = fields.index("stitle")
        qlen = fields.index("qlen")
        length = fields.index("length")
        mismatch = fields.index("mismatch")
        gap = fields.index("gaps")
        common = fields.index("common_ancestor")

        # initializes the first segment within the output
        inline = next(csvreader, None)
        current = inline[0]
        # homology score between blast and segment that takes into acount lengths of segment, length of blast hit, mismatches, and gaps
        score = (float(inline[length]) - float(inline[mismatch]) - float(inline[gap])) / float(inline[qlen])
        # stores important information within the segment
        segment.append(
            str(score) + "\t" + inline[title] + "\t" + inline[status] + "\t" + inline[source] + "\t" + inline[
                dist] + "\t" + inline[common])
        # score used for thresholding which taxa to save
        maximum = score
        inline = next(csvreader, None)

        while inline != None:
            # if contig and segment do not match previous segment
            if not inline[0] == current:
                # threshold score for keeping taxa
                thresh = maximum - args.percent
                # get the span and contig for previous contig
                contigspan = current.rsplit('_', 1)
                span = contigspan[1].split("-")
                # remove all taxa that do not meet the homology score threshold
                subset = subsetsegment(segment, thresh, graphdata, contigspan[0], span)
                # determine if segment is not a high quality match to taxa within rank of taxa of interest; if so, score segment
                keep, segscore, modsegment, justify = keepSegment(subset, args.taxaid, args.threshold, ranklin, cutoff,
                                                                  ncbi)
                # keep the segment if it is not and then save all segment information for later processing
                if keep == True:
                    vals = current
                    for hit in modsegment:
                        vals += "\t" + hit
                    data.append(vals)
                    scores.append(segscore)
                    justification.append(justify)
                    count += 1
                # reset segment to current segment and save taxa information
                segment = []
                score = (float(inline[length]) - float(inline[mismatch]) - float(inline[gap])) / float(inline[qlen])
                maximum = score
                segment.append(
                    str(score) + "\t" + inline[title] + "\t" + inline[status] + "\t" + inline[source] + "\t" + inline[
                        dist] + "\t" + inline[common])
                current = inline[0]
            else:
                # if they are the same identify, if segment meets current score threshold; save if it does
                score = (float(inline[length]) - float(inline[mismatch]) - float(inline[gap])) / float(inline[qlen])
                if score > maximum:
                    maximum = score
                if (score + args.percent) >= maximum:
                    segment.append(
                        str(score) + "\t" + inline[title] + "\t" + inline[status] + "\t" + inline[source] + "\t" +
                        inline[dist] + "\t" + inline[common])
            # read next taxa information from file
            inline = next(csvreader, None)

        # determine if necessary to keep last segment and save taxa and score information if it is
        thresh = maximum - args.percent
        contigspan = current.rsplit('_', 1)
        span = contigspan[1].split("-")
        subset = subsetsegment(segment, thresh, graphdata, contigspan[0], span)
        keep, segscore, modsegment, justify = keepSegment(subset, args.taxaid, args.threshold, ranklin, cutoff, ncbi)

        if keep == True:
            vals = current
            for hit in modsegment:
                vals += "\t" + hit
            data.append(vals)
            scores.append(segscore)
            justification.append(justify)
            count += 1

        # combine overlapping segments and then score segments
        printdata(data, virus_data, args.threshold, args.taxaid, sequencelist, coveragedata, ncbi, scores,
                  justification, properties)

        # identify if there are any segment that do not have any blast hits
        if args.unmapped:
            unmappedfile = open(args.unmapped, 'r')
            parseunmapped(unmappedfile, virus_data)
            unmappedfile.close()

        # determine if there are any deletions versus reference strain
        if args.deletion:
            deletionfile = open(args.deletion, 'r')
            parsedeletion(deletionfile, virus_data)
            deletionfile.close()

        # save combined segment information in a json format
        json.dump(virus_data, outfile, indent="\t")
        # create figures if any are desired
        if args.figname:
            createfigures(graphdata, args.figname, ncbi)
        # close files
        outfile.close()
        csvfile.close()
    else:
        # if a given rank is not found within the lineage, report to error file
        errfile.write("Could not find taxa of rank " + args.rank + " in tree of taxaid: " + args.taxaid + "\n")
    errfile.close()


def getlca(taxarank, ranklin, taxahit, ncbi):
    # identifies a pseudo least common ancestor based on cutoffs rather than absolute lca
    # lca definition: 0 - taxa within the same rank cutoff as organism of interest
    # 1- virus, 2 - nonvirus

    # default pseudo lca is non virus
    lca = 2
    # get lineage of the blast hit
    hitlin = ncbi.get_lineage(int(taxahit))
    # default absolute lca is root
    low = 1

    # if the taxa is within the rank cutoff, treat as same taxa
    if taxarank in hitlin:
        lca = 0
    # otherwise check to see if the taxa is a virus
    elif 10239 in hitlin:
        lca = 1
        
    return lca


def createfigures(graphdata, figname, ncbi):
    # creates simple charts that highlight homology score and closest taxa (by distance) to the taxa of interest
    for contig in graphdata:
        file = open(figname + "_" + contig + "_data.txt", 'w')
        start = []
        topscore = []
        toptaxa = []
        tophit = []
        for point in graphdata[contig]:
            start.append(int(point[0]))
            topscore.append(float(point[3]))
            toptaxa.append(int(point[4]))
            tophit.append(int(point[2]))
            taxaname = ncbi.get_taxid_translator([int(point[4])])
            name = "Missing_Taxa_ID"
            if int(point[4]) in taxaname:
                name = taxaname[int(point[4])]
            file.write(
                str(point[0]) + "\t" + name + "\t" + str(point[4]) + "\t" + str(point[3]) + "\t" + str(point[2]) + "\n")
        plt.plot(start, topscore, label="Top Scores")
        plt.xlabel("Genome Position")
        plt.ylabel("Sequence Score")
        plt.savefig(figname + "_" + contig + "_topscore.png")
        plt.close()
        plt.plot(start, toptaxa, label="Top Taxa")
        plt.xlabel("Genome Position")
        plt.ylabel("TaxaID")
        plt.savefig(figname + "_" + contig + "_toptaxa.png")
        plt.close()
        plt.plot(start, tophit, label="Taxa Distance")
        plt.xlabel("Genome Position")
        plt.ylabel("Distance")
        plt.savefig(figname + "_" + contig + "_tophit.png")
        plt.close()

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(start, topscore)
        ax2.plot(start, toptaxa)
        ax1.set_ylabel("Sequence Score")
        ax2.set_ylabel("TaxaID")
        ax2.set_xlabel("Genome Position")
        plt.savefig(figname + "_" + contig + "_score_taxa.png")
        plt.close()

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(start, topscore)
        ax2.plot(start, tophit)
        ax1.set_ylabel("Sequence Score")
        ax2.set_ylabel("Distance")
        ax2.set_xlabel("Genome Position")
        plt.savefig(figname + "_" + contig + "_score_hit.png")
        plt.close()
        file.close()


def subsetsegment(segment, cutoff, graphdata, contig, span):
    # removes taxa that are not within a threshold of the highest homology score

    # saves taxa information for all taxa with higher scores than threshold
    subset = []
    # arbitary maximum distance used for creating figures
    tophit = 2000000000
    # homology score and taxa of the taxa with the highest score
    topscore = 0
    toptaxa = 0

    # for each taxa in segment, determine if it is above the cutoff homology score
    # if it is save the taxa information
    for hits in segment:
        tab = hits.split("\t")
        if float(tab[0]) >= cutoff:
            subset.append(tab[0] + "\t" + tab[1] + "\t" + tab[2] + "\t" + tab[3] + "\t" + tab[4] + "\t" + tab[5])
            # when saving closest taxa for creating figures, ignore taxa without distances
            if not "CLOSEST_TAXID_NOT_FOUND" in tab[2]:
                # if distance is closer, update taxa and score
                if (int(tab[4]) < tophit) and (int(tab[4]) > -1):
                    tophit = int(tab[4])
                    topscore = float(tab[0])
                    toptaxa = tab[3]
                # if distance is the same, see if the taxa has a higher score; if so, update
                elif int(tab[4]) == tophit:
                    if float(tab[0]) > topscore:
                        topscore = float(tab[0])
                        toptaxa = tab[3]
    # save the closest taxa above threshold for figures
    if not contig in graphdata:
        graphdata[contig] = []
    graphdata[contig].append([span[0], span[1], tophit, topscore, toptaxa])

    # returns taxa information above the cutoff value for a single segment
    return subset


def printdata(data, virus_data, threshold, taxaid, sequencelist, coveragedata, ncbi, scores, justification, properties):
    # combines overlapping segments and initializes the json format
    # maximizes the score for each combined segment across the individual segments based on highest value and saves any/all unique justifications for the highest score

    # previous segment when current segment doesn't overlap
    prev = ""
    # segment or combined segment
    current = ""
    # start of segment/combined segment
    segstart = ""
    # end of segment/combined segment
    segend = ""
    # justification for the score
    justify = ""
    # score for the segment/combined segment
    segscore = 0
    # all taxa within a segment/combined segment
    total = []
    # flag for whether or not a segment has quality or coverage issues
    special = None

    # initializes the taxa to having no name and then searches and updates the names (if it has one)
    tname = "Missing_Taxa_ID"
    taxaname = ncbi.get_taxid_translator([taxaid])
    if taxaid in taxaname:
        tname = taxaname[taxaid]

    # creates the json report for target taxa of interest and run parameters
    virus_data["target_species"] = {"taxaid": taxaid, "taxa_name": tname, "rank_cutoff": properties[0],
                                    "high_quality_threshold": properties[1],
                                    "quality_difference_threshold": properties[2]}

    # for each segment within the data, combine overlapping segments and maximize the score across the combine segment
    for vals in data:
        # if first segment, initializes all parameters being saved
        if current == "":
            tab = vals.split("\t")
            name = tab[0].rsplit('_', 1)
            current = name[0]
            segment = name[1].split("-")
            segstart = int(segment[0])
            segend = int(segment[1])
            segscore = scores[0]
            justify = justification[0]
            justification.pop(0)
            scores.pop(0)
            tab.pop(0)
            total.append(tab)
        # otherwise check to see if segments overlap
        else:
            tab = vals.split("\t")
            name = tab[0].rsplit('_', 1)
            # checks to see if in the same contig
            if name[0] == current:
                segment = name[1].split("-")
                # checks to see if the segments overlap, if so combine taxa information
                if (int(segment[0]) > segstart) and (int(segment[0]) <= segend):
                    tab.pop(0)
                    total.append(tab)
                    segend = int(segment[1])
                    # maximizes score across combined segments
                    if scores[0] > segscore:
                        segscore = scores[0]
                        justify = justification[0]
                    # if multiple segment score are the max, save all unique justifications
                    elif scores[0] == segscore:
                        if not justification[0] in justify:
                            justify = justify + " and " + justification[0]
                    scores.pop(0)
                    justification.pop(0)
                else:
                    # if segments are not overlapping, format the previous segment in json format and then update all saved information for current segment
                    sequence = sequencelist[current][(int(segstart) - 1):int(segend)]
                    if coveragedata:
                        special = getspecial(current, segstart, segend, coveragedata)
                    printresults(current, segstart, segend, total, virus_data, threshold, taxaid, sequence, special,
                                 ncbi, segscore, justify)
                    prev = current
                    current = name[0]
                    segment = name[1].split("-")
                    segstart = int(segment[0])
                    segend = int(segment[1])
                    total = []
                    segscore = scores[0]
                    justify = justification[0]
                    justification.pop(0)
                    scores.pop(0)
                    tab.pop(0)
                    total.append(tab)
            else:
                # if contigs are not the same, format the previous segment in json format and then update all saved information for current segment
                sequence = sequencelist[current][(int(segstart) - 1):int(segend)]
                if coveragedata:
                    special = getspecial(current, segstart, segend, coveragedata)
                printresults(current, segstart, segend, total, virus_data, threshold, taxaid, sequence, special, ncbi,
                             segscore, justify)
                prev = current
                current = name[0]
                segment = name[1].split("-")
                segstart = int(segment[0])
                segend = int(segment[1])
                segscore = scores[0]
                egscore = scores[0]
                justify = justification[0]
                justification.pop(0)
                scores.pop(0)
                total = []
                tab.pop(0)
                total.append(tab)

    # format final segment in json format
    if data:
        sequence = sequencelist[current][(int(segstart) - 1):int(segend)]
        if coveragedata:
            special = getspecial(current, segstart, segend, coveragedata)
        printresults(current, segstart, segend, total, virus_data, threshold, taxaid, sequence, special, ncbi, segscore,
                    justify)


def printresults(current, segstart, segend, total, virus_data, threshold, taxaid, sequence, special, ncbi, segscore,
                 justify):
    # formats individual or combined segments into json format by identifying and categorizing different properties of the segment
    # breaks up groups into low homology within taxa rank, high/low homology outside rank, blacklist, taxaid not found

    # segment counts containing low homology taxa within rank cutoff
    lowsame = {}
    # taxa information for high homology nonblacklist taxa that are outside the rank cutoff
    highdif = {}
    # taxa information for low homology nonblacklist taxa that are outside the rank cutoff
    lowdif = {}
    # taxa information for blacklist taxa
    blacklist = {}
    # taxa information for taxa without ids
    other = {}
    # number of individual segments within the combined segment (can just be an individual segment)
    segnum = 0
    # number of segments that contain low homology within rank
    ls = 0
    # number of segments that contain taxa without ids
    ot = 0
    # number of segments that contain high homology outside rank (nonblacklist)
    hd = 0
    # number of segments that contain low homology outside rank (nonblacklist)
    ld = 0
    # number of segments that contain blacklist
    bl = 0
    # booleans for each above segment count used to determine if an individual segment has a taxa with that feature
    bot = False
    bls = False
    bhd = False
    bld = False
    bbl = False

    for val in total:
        # initializes each individual segment to the first taxa and having no features
        num = 0
        bls = False
        bot = False
        bhd = False
        bld = False
        bbl = False
        # while there are still taxa in the individual segment
        while num < len(val):
            # if it doesn't have a taxa id, save to the other array
            if "CLOSEST_TAXID_NOT_FOUND" in val[2 + num]:
                # if taxa has not been previously found, add to array: number of segments the taxa is in, name of taxa, pseudo lca, last segment in, homology score, each segement present
                if not val[3 + num] in other:
                    other[val[3 + num]] = [1, val[1 + num], val[4 + num], segnum, val[0 + num], str(segnum + 1)]
                    bot = True
                else:
                    arr = other[val[3 + num]]
                    # if homology score is higher than current highest, save that score
                    if float(val[0 + num]) > float(arr[4]):
                        arr[4] = val[0 + num]
                    # otherwise check to see if taxa names match, append if they do not
                    if not val[1 + num] in arr[1]:
                        arr[1] = arr[1] + ";" + val[1 + num]
                    # check to see if not seen in this segment before; if so, increase number of segments the taxa is in
                    if segnum > arr[3]:
                        arr[3] = segnum
                        arr[0] += 1
                        arr[5] += "; " + str(segnum + 1)
                        bot = True
                    other[val[3 + num]] = arr
            # if its blacklist, save information to blacklist array
            if "BLACKLISTED" in val[2 + num]:
                # if taxa has not been previously found, add to array: number of segments the taxa is in, name of taxa, pseudo lca, last segment in, homology score, each segement present
                if not val[3 + num] in blacklist:
                    blacklist[val[3 + num]] = [1, val[1 + num], val[4 + num], segnum, val[0 + num], str(segnum + 1)]
                    bbl = True
                else:
                    arr = blacklist[val[3 + num]]
                    # if homology score is higher than current highest, save that score
                    if float(val[0 + num]) > float(arr[4]):
                        arr[4] = val[0 + num]
                    # otherwise check to see if taxa names match, append if they do not
                    if not val[1 + num] in arr[1]:
                        arr[1] = arr[1] + ";" + val[1 + num]
                    # check to see if not seen in this segment before; if so, increase number of segments the taxa is in
                    if segnum > arr[3]:
                        arr[3] = segnum
                        arr[0] += 1
                        arr[5] += "; " + str(segnum + 1)
                        bbl = True
                    blacklist[val[3 + num]] = arr
            else:
                # check to see if taxa is within rank cutoff
                if (val[4 + num] == "0"):
                    # if taxa has not been previously found, add to array: number of segments the taxa is in, name of taxa, pseudo lca, last segment in, homology score, each segement present
                    if not val[3 + num] in lowsame:
                        lowsame[val[3 + num]] = [1, val[1 + num], val[4 + num], segnum, val[0 + num], str(segnum + 1)]
                        bls = True
                    else:
                        arr = lowsame[val[3 + num]]
                        # if homology score is higher than current highest, save that score
                        if float(val[0 + num]) > float(arr[4]):
                            arr[4] = val[0 + num]
                        # otherwise check to see if taxa names match, append if they do not
                        if not val[1 + num] in arr[1]:
                            arr[1] = arr[1] + ";" + val[1 + num]
                        # check to see if not seen in this segment before; if so, increase number of segments the taxa is in
                        if segnum > arr[3]:
                            arr[3] = segnum
                            arr[0] += 1
                            arr[5] += "; " + str(segnum + 1)
                            bls = True
                        lowsame[val[3 + num]] = arr
                else:
                    # otherwise check to see if the different taxa has a low or high quality homology score
                    if float(val[0 + num]) <= threshold:
                        # if taxa has not been previously found, add to array: number of segments the taxa is in, name of taxa, pseudo lca, last segment in, homology score, each segement present
                        if not val[3 + num] in lowdif:
                            lowdif[val[3 + num]] = [1, val[1 + num], val[4 + num], segnum, val[0 + num],
                                                    str(segnum + 1)]
                            bld = True
                        else:
                            arr = lowdif[val[3 + num]]
                            # if homology score is higher than current highest, save that score
                            if float(val[0 + num]) > float(arr[4]):
                                arr[4] = val[0 + num]
                            # otherwise check to see if taxa names match, append if they do not
                            if not val[1 + num] in arr[1]:
                                arr[1] = arr[1] + "\t" + val[1 + num]
                            # check to see if not seen in this segment before; if so, increase number of segments the taxa is in
                            if segnum > arr[3]:
                                arr[3] = segnum
                                arr[0] += 1
                                arr[5] += "; " + str(segnum + 1)
                                bld = True
                            lowdif[val[3 + num]] = arr
                    else:
                        # if taxa has not been previously found, add to array: number of segments the taxa is in, name of taxa, pseudo lca, last segment in, homology score, each segement present
                        if not val[3 + num] in highdif:
                            highdif[val[3 + num]] = [1, val[1 + num], val[4 + num], segnum, val[0 + num],
                                                     str(segnum + 1)]
                            bhd = True
                        else:
                            arr = highdif[val[3 + num]]
                            # if homology score is higher than current highest, save that score
                            if float(val[0 + num]) > float(arr[4]):
                                arr[4] = val[0 + num]
                            # otherwise check to see if taxa names match, append if they do not
                            if not val[1 + num] in arr[1]:
                                arr[1] = arr[1] + "\t" + val[1 + num]
                            # check to see if not seen in this segment before; if so, increase number of segments the taxa is in
                            if segnum > arr[3]:
                                arr[3] = segnum
                                arr[0] += 1
                                arr[5] += "; " + str(segnum + 1)
                                bhd = True
                            highdif[val[3 + num]] = arr
            # update the position in the array to the next taxa entry for the segment
            num += 6

        # if a given category is seen in the current segment, increase the count total for the category
        if bot is True:
            ot += 1
        if bbl is True:
            bl += 1
        if bls is True:
            ls += 1
        if bld is True:
            ld += 1
        if bhd is True:
            hd += 1
        # set the current segment number to the next segment in the combined segment
        segnum += 1

    # initializes the segment to contain start, end, number of segments, nucleic acid sequence, maximized score, and justification
    posinfo = {}
    posinfo["start"] = segstart
    posinfo["end"] = segend
    posinfo["number_segments"] = segnum
    posinfo["sequence"] = sequence
    posinfo["score"] = segscore
    posinfo["justification"] = justify

    # determines if the combined segment contains low coverage basepairs
    lowcov = False
    specialstring = ""
    if special:
        posinfo["special"] = []
        num = 1
        for item in special:
            sets = item.split(",")
            if "LowCoverage" in sets[3]:
                lowcov = True
            posinfo["special"].append({"type": sets[3], "start": int(sets[0]), "end": int(sets[1])})
            num += 1

    # if there are any low homology hits within the taxa rank of interest, create json entry
    if ls > 0:
        posinfo["low_homology_within_target_rank"] = {"number_segments": ls}

    # if there are any high quality nonblacklist hits to taxa outside that rank, create json entry
    if hd > 0:
        taxainfo = printarr(highdif, ncbi)
        if len(highdif) > 10:
            posinfo["high_quality_different_taxa"] = {"number_taxa": len(highdif),
                                                      "comment": "showing most prevalent taxa", "taxa": taxainfo}
        else:
            posinfo["high_quality_different_taxa"] = {"number_taxa": len(highdif), "comment": "showing all taxa",
                                                      "taxa": taxainfo}

    # if there are any low quality nonblacklist hits to taxa outside that rank, create json entry
    if ld > 0:
        taxainfo = printarr(lowdif, ncbi)
        if len(lowdif) > 10:
            posinfo["low_quality_different_taxa"] = {"number_taxa": len(lowdif),
                                                     "comment": "showing most prevalent taxa", "taxa": taxainfo}
        else:
            posinfo["low_quality_different_taxa"] = {"number_taxa": len(lowdif), "comment": "showing all taxa",
                                                     "taxa": taxainfo}

    # if there are any blacklist hits, create json entry
    if bl > 0:
        taxainfo = printarr(blacklist, ncbi)
        if len(blacklist) > 10:
            posinfo["engineered_or_lab_taxa"] = {"number_taxa": len(blacklist), "comment": "showing most prevalent taxa",
                                         "taxa": taxainfo}
        else:
            posinfo["engineered_or_lab_taxa"] = {"number_taxa": len(blacklist), "comment": "showing all taxa", "taxa": taxainfo}

    # if there is any low coverage regions in the combined segment, update the score, justification, and position in the json
    if lowcov:
        if not "low_coverage" in virus_data:
            virus_data["low_coverage"] = {}
        if not current in virus_data["low_coverage"]:
            virus_data["low_coverage"][current] = []
        posinfo["score"] = 1
        posinfo["justification"] = "low coverage region"
        virus_data["low_coverage"][current].append(posinfo)
    else:
        # otherwise the combined segment is a positive hit
        if not "regions_of_interest" in virus_data:
            virus_data["regions_of_interest"] = {}
        if not current in virus_data["regions_of_interest"]:
            virus_data["regions_of_interest"][current] = []
        virus_data["regions_of_interest"][current].append(posinfo)


def printarr(sets, ncbi):
    # prints taxaid counts distance taxanames
    taxainfo = []
    # if there are 10 or less taxa associated with each taxa category, print all taxa
    if len(sets) <= 10:
        for key in sets:
            arr = sets[key]
            taxaname = ncbi.get_taxid_translator([int(key)])
            name = "Missing_Taxa_ID"
            taxatype = "NA"
            if arr[2] == "2":
                taxatype = "nonvirus"
            elif arr[2] == "1":
                taxatype = "virus"
            if int(key) in taxaname:
                name = taxaname[int(key)]
            taxainfo.append({"taxaid": int(key), "name": name, "segments": arr[0], "present_in_segment(s)": arr[5],
                             "max_homology_score": arr[4], "taxa_type": taxatype})
    else:
        # otherwise find the highest homology score
        maximum = 0
        taxa = []
        for key in sets:
            arr = sets[key]
            if int(arr[0]) > maximum:
                maximum = int(arr[0])
        for key in sets:
            arr = sets[key]
            if int(arr[0]) == maximum:
                taxa.append([int(arr[2]), arr[1], key, arr[4], arr[2], arr[5]])
        sortedtaxa = sorted(taxa, key=lambda x: x[0])

        tot = 0
        # of the taxa with the highest homology score, print at most 10
        for key in sortedtaxa:
            if tot < 10:
                taxaname = ncbi.get_taxid_translator([int(key[2])])
                name = "Missing_Taxa_ID"
                taxatype = "NA"
                if key[4] == "2":
                    taxatype = "nonvirus"
                elif key[4] == "1":
                    taxatype = "virus"
                if int(key[2]) in taxaname:
                    name = taxaname[int(key[2])]
                taxainfo.append(
                    {"taxaid": int(key[2]), "name": name, "segments": maximum, "present_in_segment(s)": key[5],
                     "max_homology_score": key[3], "taxa_type": taxatype})
                tot += 1
            else:
                break

    return taxainfo


def keepSegment(segment, taxaid, threshold, ranklin, taxarank, ncbi):
    # determines if the segment is not within the rank cutoff; if it is not, score and save the segment

    # defaults to keeping the segment
    keep = True
    # if there is a high quality blacklist taxa within the segment
    HBL = False
    # if there is a low quality blacklist taxa within the segment
    LBL = False
    # if there is a high quality virus taxa within the segment
    hvirus = False
    # if there is a low quality virus taxa within the segment
    lvirus = False
    # if there is a high quality nonvirus taxa within the segment
    hnv = False
    # if there is a low quality nonvirus taxa within the segment
    lnv = False
    # if there is a low quality hit within the taxa rank cutoff within the segment
    lsv = False
    # if there is a high quality taxa without taxa id within the segment
    hcnf = False
    # if there is a low quality taxa without taxa id within the segment
    lcnf = False
    # default to region of loq confidence
    score = 4
    # reason for the score
    justify = ""
    # saves taxa information replacing distance with a pseudo lca
    modsegment = []

    # for each taxa, determine where it falls within the taxa types of interest and whether or not it is a high quality hit
    for item in segment:
        name = item.split("\t")
        # defaults to within taxa rank
        lca = 0
        if "CLOSEST_TAXID_NOT_FOUND" in name[2]:
            # taxa without id has no lca
            lca = 3
            if float(name[0]) > threshold:
                hcnf = True
            else:
                lcnf = True
        elif (int(name[3]) == taxaid):
            # determines if the taxa id is the same as the taxa of interest; if it is a high quality hit, do not keep segment
            if float(name[0]) > threshold:
                keep = False
                break;
            else:
                lsv = True
        elif "BLACKLISTED" in name[2]:
            lca = getlca(taxarank, ranklin, name[3], ncbi)
            if float(name[0]) > threshold:
                HBL = True
            else:
                LBL = True
        else:
            # gets a pseudo lca
            lca = getlca(taxarank, ranklin, name[3], ncbi)
            # determines if the taxa is within the rank cutoff; if it is a high quality hit, do not keep segment
            if lca == 0:
                if float(name[0]) > threshold:
                    keep = False
                    break;
                else:
                    lsv = True
            # determines if the taxa is a virus
            elif lca == 1:
                if float(name[0]) > threshold:
                    hvirus = True
                else:
                    lvirus = True
            # otherwise it is a nonvirus
            else:
                if float(name[0]) > threshold:
                    hnv = True
                else:
                    lnv = True
        # saves taxa information replacing distance with pseudo lca
        modsegment.append(name[0] + "\t" + name[1] + "\t" + name[2] + "\t" + name[3] + "\t" + str(lca) + "\t" + name[5])

    # minimizes individual segment score (ie. always chooses value closest to a if it exists)
    # first based on: a) within taxa rank, b) virus, c) nonvirus, d) blacklist, e) no taxa id
    # after based on: a) high quality, b) low quality
    # adds a statement justify why the score was chosen
    if keep == True:
        if lsv:
            score = 2
            justify = "low homology within taxa rank"
        elif hvirus:
            score = 5
            justify = "high quality hit to virus"
        elif lvirus:
            score = 2
            justify = "low quality hit to virus"
        elif hnv:
            score = 9
            justify = "high quality hit to nonvirus"
        elif lnv:
            score = 4
            justify = "low quality hit to nonvirus"
        elif HBL:
            score = 10
            justify = "high quality hit to engineered/lab taxa"
        elif LBL:
            score = 5
            justify = "low quality hit to engineered/lab taxa"
        elif hcnf:
            score = 9
            justify = "high quality hit to organism without taxa"
        elif lcnf:
            score = 4
            justify = "low quality hit to organism without taxa"

    return keep, score, modsegment, justify


def getspecial(current, segstart, segend, coveragedata):
    # identifies if a given segment has any coverage or quality issues
    special = []
    sets = None
    if current in coveragedata:
        sets = coveragedata[current]

    start = int(segstart)
    end = int(segend)

    if sets:
        # for all areas within coverage/quality issues, determine if they exist in the current segment
        for item in sets:
            elem = item.split(",")
            estart = int(elem[0])
            eend = int(elem[1])

            if (estart >= start) and (estart <= end):
                special.append(item)
            elif ((eend >= start) and (eend <= end)):
                special.append(item)
            elif ((estart < start) and (eend > end)):
                special.append(item)

    return special


def parsedeletion(deletionfile, virus_data):
    # parses and scores regions that exist within the reference but do not exist within the assembly
    inline = deletionfile.readline()
    count = 1
    while inline:
        inline = inline.strip()
        if not "#" in inline:
            span = inline.split(" | ")
            if int(span[4]) > 0:
                if not "deletion" in virus_data:
                    virus_data["deletion"] = {}
                name = span[0]
                if not name in virus_data["deletion"]:
                    virus_data["deletion"][name] = []
                delinfo = {}
                delinfo["assembly_affected"] = span[1]
                delinfo["start"] = span[2]
                delinfo["end"] = span[3]
                delinfo["length"] = span[4]
                delinfo["score"] = 2
                delinfo["justification"] = "deletion"
                virus_data["deletion"][name].append(delinfo)
                count += 1
        inline = deletionfile.readline()


def parseunmapped(unmappedfile, virus_data):
    # parses and scores segments within that assembly that are not low complexity and do not have any blast hits
    inline = unmappedfile.readline()

    while inline:
        inline = inline.strip()
        if ">" in inline:
            if not "novel" in virus_data:
                virus_data["novel"] = {}
            current = inline.replace(">", "")
            contigspan = current.rsplit('_', 1)
            span = contigspan[1].split("-")
            name = contigspan[0]
            if not name in virus_data["novel"]:
                virus_data["novel"][name] = []
            uninfo = {}
            uninfo["start"] = span[0]
            uninfo["end"] = span[1]
            inline = unmappedfile.readline()
            inline = inline.strip()
            uninfo["sequence"] = inline
            uninfo["score"] = 1
            uninfo["justification"] = "no hits to known taxa"
            virus_data["novel"][name].append(uninfo)
        inline = unmappedfile.readline()


def parsecoverage(coveragefile):
    # identifies and saves assembly positions that contain any coverage or quality issues
    coveragedata = {}

    line = coveragefile.readline()
    name = ""

    while line:
        line = line.strip()
        sets = line.split(";")
        count = 0
        features = []
        for item in sets:
            if count == 0:
                name = item.split(",")
                count += 1
            else:
                features.append(item)
        coveragedata[name[0]] = features
        line = coveragefile.readline()

    return coveragedata


def parsefasta(fastafile):
    # saves sequence for each contig used in the analysis
    sequencelist = {}
    sequence = ""
    name = ""

    line = fastafile.readline()

    while line:
        line = line.strip()
        if ">" in line:
            if name == "":
                name = line.replace(">", "")
            else:
                sequencelist[name] = sequence
                name = line.replace(">", "")
                sequence = ""
        else:
            sequence += line
        line = fastafile.readline()

    sequencelist[name] = sequence

    return sequencelist


if __name__ == "__main__":
    sys.exit(main())

