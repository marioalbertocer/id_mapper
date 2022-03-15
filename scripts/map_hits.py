import os, sys
import pandas as pd

"""
clustering() is a function that generates a cluster of blast hits based on the contig id (node), the strand (positive 
or negative), and the element (chromosome, psymA or psymB). For instance, a cluster of hits of the positive strand of
the contig/node 20 against chromosomic sequences of E. memiloti. Inputs for this function are the info from the blast
report (f), the node or contig id (node), the first nucleotide position of the hit (smin), the strand, and the element.  
"""

def clustering(f, node, smin, strand, element):
	cluster = []
	for line in f:
		line = line.strip()
		nodeIN = line.split(",")[1].split("_length_")[0] 
		if (line.split(",")[0]).split("_")[0] == element:
			if nodeIN == node:
				al_length = int(line.split(",")[3])
				if (al_length >= 50): 
					if line.endswith(strand):
						min = int(line.split(",")[smin])
						if len(cluster) == 0:
							cluster.append(line)
						else:
							ind = 0
							n = "c"
							for i in cluster:
								ind += 1
								mmin = int(i.split(",")[smin])
								if min <= mmin:
									cluster.insert((ind-1) , line)
									n = "nc"
									break
							if n == "c":
								cluster.append(line)
	return cluster

"""
add_hits() is a function that calculates the total coverage of hits in a cluster of hits. It considers overlapping of
hits positions. For instance, if the fist, second, and third hit in a cluster are between the positions 1-15, 10-20, y 
30-40, respectively. Then, the total coverage of these hits  should be 30 (15 from the first hit plus 5 for the 
non-overlapping portion of the second hit and the 10 of the third hit). The inputs for this function are a cluster of 
hits and the number of the columns that contain the fist and the last postions.   
"""

def add_hits(cluster, smin, smax):
	positions_MnMx = {}
	counter = 0
	for line in cluster:
		positions_MnMx[counter] = [int(line.split(",")[smin]), int(line.split(",")[smax])]
		counter += 1
	df_positions = pd.DataFrame(positions_MnMx) # convert positions_MnMx dictionary to a dataframe
	dft_positions = df_positions.T # swap columns and rows of the dataframe -- now index is raws and smin/smax is columns. 
	dft_positions = dft_positions.sort_values(by = [0,1], ascending = [True, True], na_position = 'first') # sorts both columns 
	dfl_positions = dft_positions.values.tolist() # converts the whole info into a list like "[[1, 2], [3, 4], [5, 6]]"
	total = 0
	mn = 0
	mx = 0
	new_mn = 0
	new_mx = 0
	for i in dfl_positions:
		new_mn = i[0]
		new_mx = i[1]
		if total != 0:
			if new_mx > mx:
				if new_mn < mx:
					mn = mx
					mx = new_mx
				else:
					mn = new_mn			
					mx = new_mx
			else:
				continue
		else:
			mn = new_mn			
			mx = new_mx
		total += mx - mn
	return total

"""
classify()

"""

def classify(coverage_e_s, l_cutoff, d_ratio):
	sortedCounts = sorted(coverage_e_s.items(), key=lambda x: x[1], reverse=True) # sorts (descending order) the library by hit coverage length
	bestID = ""
	div1n2 = ""
	lengthCut = ""
	
	if (sortedCounts[0])[1] > l_cutoff: # if longest coverage is higher than 3000
		lengthCut = "Y"  # it passes the length coverage criterion  
	else:
		lengthCut = "N"
	
	if (sortedCounts[1])[1] != 0: # checking on the next longest hit coverage
		if ((sortedCounts[1])[0]).split("_")[0] != ((sortedCounts[0])[0]).split("_")[0]:
			div1n2 =  str(round(((sortedCounts[0])[1] / (sortedCounts[1])[1]), 2)) # comparing the two values (two longest coverages)
			if float(div1n2) >= d_ratio:
				bestID = (sortedCounts[0])[0] # If the longest is more than twice a long as the second longest
			else:
				bestID = "undetermined"
		else:
			div1n2 = "same"  # if the two longest coverages are from the same element (different strand)
			bestID = (sortedCounts[0])[0]
	else:
		div1n2 = "NA" # if there is NOT a next longest (hits only from a strand of an element)
		bestID = (sortedCounts[0])[0]
	return bestID, div1n2, lengthCut