from Bio import SeqIO
import os,sys 

def sortSeqs(elements, path2out, assemblies, rep):
	out_categories = elements + ["undetermined", "too_small"]
	
	#	paths_to_sorted_sequences
	paths = []
	paths.append(path2out + "sorted/")
	for element in out_categories: paths.append(path2out + "sorted/" + element)
	for path_i in paths:
		if not os.path.exists(path_i): os.system("mkdir " + path_i)
	
	processed = 0
	for al in assemblies:
		file_name, strain = map(lambda x: x.strip(), al.split(","))
		sequences = {}
		names = {}
		e_fas = {}
		for element in out_categories: e_fas[element] = []
		for rec in SeqIO.parse(file_name, "fasta"):
			nodeID = (rec.id).split("_length")[0]
			names[nodeID] = str(rec.id)
			sequences[nodeID] = str(rec.seq)
		for rl in rep:
			rl = rl.strip()
			if rl.split(",")[0] == strain:
				idn = (rl.split(",")[-3]).split("_")[0]
				node = rl.split(",")[1]
				if rl.endswith("Y"):
					if idn in elements + ["undetermined"]: e_fas[idn].append([names[node], sequences[node]])
				else:
					e_fas["too_small"].append([names[node], sequences[node]])
		for element in out_categories:
			sorted_seqs = open(paths[0] + element + "/" + file_name.split("/")[-1], "w")
			for e, name_seqs in e_fas.items():
				if element == e:
					if name_seqs != []:
						for name_seq in name_seqs:
							sorted_seqs.write(">%s\n%s\n" % (name_seq[0], name_seq[1]))
			sorted_seqs.close()
		processed += 1
	finished = "n"
	if processed == len(assemblies): finished = "y"
	return finished
	
	