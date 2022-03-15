from Bio import SeqIO
import os,sys 

def sortSeqs(elements, path2out, mapping, rep):
	out_categories = elements + ["undetermined", "too_small"]

	#	paths_to_sorted_sequences
	paths = []
	paths.append(path2out + "sorted/")
	for element in out_categories: paths.append(path2out + "sorted/" + element)
	for path_i in paths:
		if not os.path.exists(path_i): os.system("mkdir " + path_i)
	
	processed = 0
	for ml in mapping:
		if ml.startswith("assembly"):
			file_name, path2file = map(lambda x: x.strip(), ml.split(",")[1:3])
			sequences = {}
			names = {}
			e_fas = {}
			for element in out_categories: e_fas[element] = []
			for rec in SeqIO.parse(path2file + "/" + file_name, "fasta"):
				nodeID = (rec.id).split("_length")[0]
				names[nodeID] = str(rec.id)
				sequences[nodeID] = str(rec.seq)
			for rl in rep:
				rl = rl.strip()
				if rl.split(",")[0] == file_name.split("_")[0]:
					idn = (rl.split(",")[9]).split("_")[0]
					node = rl.split(",")[1]
					if rl.endswith("Y"):
						if idn in elements + ["undetermined"]: e_fas[idn].append([names[node], sequences[node]])
					else:
						e_fas["too_small"].append([names[node], sequences[node]])
			for element in out_categories:
				sorted_seqs = open(paths[0] + element + "/" + file_name, "w")
				for e, name_seqs in e_fas.items():
					if element == e:					
						if name_seqs != []:
							for name_seq in name_seqs:
								sorted_seqs.write(">%s\n%s\n" % (name_seq[0], name_seq[1]))
				sorted_seqs.close()
			processed += 1
	finished = "n"
	if processed == len(list(filter(lambda score: score.startswith("assembly"), mapping))): finished = "y"
	return finished
	
	