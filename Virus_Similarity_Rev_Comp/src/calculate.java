import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

/* A program designed to take as input a subset of the virus pairs with non-zero interactions.
 * Kmers for each virus are precreated so that they do not have to be regenerated for each comparison.
 * A quantitative proportional similarity on kmer counts is used to assess the similarity between two viruses.
 */

public class calculate {

	public static void main(String[] args) throws IOException {
		// file containing a unique subset of viruses to calculate similarity scores
		BufferedReader unique = new BufferedReader(new FileReader(args[0]));
		// file containing all kmers for the viruses of interest
		BufferedReader kmerFile = new BufferedReader(new FileReader(args[1]));
		// file containing the calculated similarity scores for the subset of viruses
		BufferedWriter outfile = new BufferedWriter(new FileWriter(args[2]));
		// hash table that contains all of the kmers for a given virus
		HashMap<String,String[]> kmerTable = new HashMap<String,String[]>(10000000);
		// class that identifies similarity scores between two viruses
		HashSimilarity similarity;
		// string that reads file lines and is initialized to contain the first line of the kmer file
		String inline = kmerFile.readLine();
		// array that contains the initial delimited information for each line
		String[] intab;
		// array that contains all the kmers for a given virus
		String[] kmers;
		// set of kmers for one of the viruses in the similarity calculation
		String[] interaction;
		// set of kmers for the other viruses in the similarity calculation
		String[] interaction2;
		// the resulting similarity score
		double score;
		
		// reads in the kmer file and stores all of the virus-kmer pairs in a hash table
		// each line represents a virus and all of its kmers
		while (inline != null)
		{
			if (!inline.equals(""))
			{
				// delimiter that separates the virus from the kmers
				intab = inline.split("%%");
				// delimitter that separates the kmers from each other
				kmers = intab[1].split("\t");
				// entry where the key is a virus and the object is a set of kmers
				kmerTable.put(intab[0], kmers);
			}
			inline = kmerFile.readLine();
		}
		
		// reads in the first line of the file containing unique virus interactions
		inline = unique.readLine();
		
		// identifies similarity scores for each virus pair
		while (inline != null)
		{
			if (!inline.equals(""))
			{
				// splits the two viruses
				intab = inline.split("\t");
				
				// creates a new object for identifying the similarity scores for each virus pair 
				similarity = new HashSimilarity();
				
				// gets the kmer profiles for the both of the virus
				interaction = kmerTable.get(intab[0]);
				interaction2 = kmerTable.get(intab[1]);

				// generates the kmer profile for the first virus
				for (int j = 0; j < interaction.length; j++) {
					similarity.add(0, interaction[j]);
				}

				// generates the kmer profile for the second virus
				for (int j = 0; j < interaction2.length; j++) {
					similarity.add(1, interaction2[j]);
				}

				// calculates the similarity score between the two viruses
				score = similarity.quantSimilarity();
				
				// if the similarity score isn't null, then save the score to the file
				if (score > -1.0)
				{
					outfile.write(intab[0] + "\t" + intab[1] + "\t" + score + "\n");
				}

				// deletes the object containing the kmer profiles for the two viruses
				similarity = null;
			}
			
			inline = unique.readLine();
		}
		
		// closes all of the files
		outfile.close();
		unique.close();
		kmerFile.close();

	}

}
