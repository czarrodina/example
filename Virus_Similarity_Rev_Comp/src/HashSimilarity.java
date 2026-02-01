import java.util.HashMap;

/* An object that contains the kmer profiles for both of the virus that are being compared.
 * It also calculates a quantitative proportional similarity based up the kmer counts for each virus
 * The kmer count contains the number of instances unique kmers show up in either virus
 */

public class HashSimilarity {
	// hash table containing all the unique kmers and their counts for both viruses
	HashMap<String,Integer[]> counts;
	
	// initializes the object to contain a hash table of 10 million keys
	public HashSimilarity()
	{
		counts = new HashMap<String,Integer[]>(10000000);
	}
	
	// adds kmers to the hash table where the number is the virus (either 0 or 1) and the kmer is 
	// the string of nucleotides associated with the virus
	public void add(int number, String kmer)
	{
		// creates the reverse complement for the kmer
		String rev = revComp(kmer);
		// a integer containing the counts for both viruses for the kmer of interest
		Integer[] vals;
		
		// if the kmer alread exists within the hash table then update the count values
		if (counts.containsKey(kmer))
		{
			// initializes the array and then gets the previous count values
			vals = new Integer[2];
			vals = counts.get(kmer);
			
			// if its the first virus, update the first virus' count value
			if (number == 0)
			{
				vals[0]++;
			}
			
			// other wise update the secod virus' count value
			else
			{
				vals[1]++;
			}
			
			// updates the hash table to reflect the updated counts
			counts.replace(kmer, vals);
		}
		
		// if the kmer was not in the hash table, see if the reverse complement has already exists
		else if (counts.containsKey(rev))
		{
			// initializes the array and then gets the previous count values
			vals = new Integer[2];
			vals = counts.get(rev);
			
			// if its the first virus, update the first virus' count value
			if (number == 0)
			{
				vals[0]++;
			}
			
			// other wise update the secod virus' count value
			else
			{
				vals[1]++;
			}
			
			// updates the hash table to reflect the updated counts
			counts.replace(rev, vals);
		}
		
		// if the kmer or its reverse complement is not in the hash table, add the kmer to the table
		else
		{
			// if its the first virus
			if (number ==0)
			{
				// initializes the kmer to have a one count for the first viruses and a zero count for the second
				vals = new Integer[2];
				vals[0] = 1;
				vals[1] = 0;
				
				// adds the kmer to the hash table
				counts.put(kmer, vals);
			}
			
			// otherwise if its the second virus
			else
			{
				// initializes the kmer to have a zero count for the first viruses and a one count for the second
				vals = new Integer[2];
				vals[0] = 0;
				vals[1] = 1;
				
				// adds the kmer to the hash table
				counts.put(kmer, vals);
			}
		}
	}
	
	// generates the reverse complement for a given kmer
	private String revComp(String forward)
	{
		// initializes the reverse complement to be empty
		String rc = "";
		
		// start at the end of the kmer to generate the reverse complement
		for (int i = forward.length() -1; i > -1; i--)
		{
			// gets the ith nucleotide and then updates to its reverse
			char cur = forward.charAt(i);
			if (cur == 'a')
			{
				rc += "t";
			}
			else if (cur == 'c')
			{
				rc += "g";
			}
			else if (cur == 't')
			{
				rc += "a";
			}
			else if (cur == 'g')
			{
				rc += "c";
			}
		}
		
		// returns the reverse complement
		return rc;
	}
	
	// calculates the quantitative similarity of kmer profiles between two viruses
	public double quantSimilarity()
	{
		// similarity value
		double sim = 0.0;
		// array containing the differences in kmer counts and the total kmer counts for both viruses
		int[] value = new int[2];
		
		// initializes the numerator to 0
		value[0] = 0;
		// initializes the denominator to 0
		value[1] = 0;
		// counts for a given kmer for both virus
		Integer[] val;
		// alterations to the numerator and denominator from the dissimilarity calculation
		int[] add;
		
		// for each unique kmer
		for (String key : counts.keySet())
		{
			// get the counts for both viruses
			val = counts.get(key);
			// calculate the contribution to the similarity scores for the kmer
			add = dissimilarity(val);
			
			// update the numerator and denominator of the similarity score
			value[0] += add[0];
			value[1] += add[1];
		}
		
		// if there is no kmer overlap between the viruses, return an error flag
		if (value[0] == value[1])
		{
			return -10.0;
		}
		
		// calculate the similarity score from the contributions of the dissimilarity scire
		sim = 1 - (value[0] / (value[1] * 1.0));
		
		// returns the resulting similarity score
		return sim;
	}
	
	// calculates the contribution of the kmer counts to the dissimilarity score
	private int[] dissimilarity(Integer[] values)
	{
		// contribution to the numerator and denominator for the dissimilarity score
		int[] dissim = new int[2];
		
		// differences in the kmer counts for both viruses
		dissim[0] = Math.abs(values[0] - values[1]);
		// sum of the kmer counts for both viruses
		dissim[1] = values[0] + values[1];
		
		// returns the contribution of the kmer counts to the dissimilarity score
		return dissim;
	}
}
