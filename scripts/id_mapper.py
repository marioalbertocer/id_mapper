import os,sys
from blast_step import *
from map_hits import *
from sorting_step import sortSeqs

"""
The function getParameters() get parameter for running the whole process through the parameters file that the user 
submits when calling "id_mapper.py parameters file mapping_file"
"""

def getParameters(parFile):
	parameters = []
	for lin in parFile:
		if lin.startswith("evalue"): parameters.append(lin.split(":")[1].strip())
		if lin.startswith("percent of identity"): parameters.append(lin.split(":")[1].strip())
		if lin.startswith("number of threads"): parameters.append(int(lin.split(":")[1].strip()))
		if lin.startswith("output path"): parameters.append(lin.split(":")[1].strip())		
		if lin.startswith("elements"): parameters.append(list(map(lambda x: x.strip(), (lin.split(":")[1].split(",")))))
		if lin.startswith("length cutoff"): parameters.append(int(lin.split(":")[1].strip()))
		if lin.startswith("decision ratio"): parameters.append(float(lin.split(":")[1].strip()))
	return parameters

def getLengths(f): # takes the information of a Blast report (f)
	nodeLengths = {}
	for fLine in f:
		fLine = fLine.strip()
		if "," in fLine:
			node = (fLine.split("_length_")[0]).split(",")[1]
			length = (fLine.split("_length_")[1]).split("_")[0]
			nodeLengths[node] = length
	return nodeLengths
		
def main():
	mode = 'n'
	if len(sys.argv) >= 4: mode = sys.argv[3]
	print("\n(1) setting variables and creating foldes and files")
	parameters_file = open(sys.argv[1], "r").readlines()
	mapping_file = open(sys.argv[2], "r").readlines()
	parameters = getParameters(parameters_file)
	output_dir = parameters[3] + "/"
	elements = parameters[4]
	os.system("mkdir " + output_dir)
	rep_out = open(output_dir + "id_report.csv", "w")
	strands_positions = {}
	strands_positions["minus"] = [9, 8]
	strands_positions["plus"] = [8, 9]
	for l in parameters_file: print(l.strip())
	print("\n(1) Done")

	print("\n(2) running the Blast step")
	set_Blast_inputs = setBlastInputs(mapping_file, elements, output_dir)
	concaFile = set_Blast_inputs[0]
	assemblies = set_Blast_inputs[1]
	blast_reports = runBlast(parameters, concaFile, assemblies)
	print("\n(2) Done")
	
	print("\n(3) mapping Blast hits")
	for repFile in blast_reports:
		rep = open(repFile.split(",")[0], "r").readlines()
		sps = repFile.split(",")[1]
		nodeLengths = getLengths(rep)
		coverage_e_s = {}
		for node, nodeLength in nodeLengths.items():
			rep_node = ("%s,%s,%s" % (sps,node, nodeLength))
			for element in elements:
				for strand, positions in strands_positions.items():
					cluster_i = clustering(rep, node, positions[0], strand, element) 
					if cluster_i == []:
						total = 0
					else:
						total = add_hits(cluster_i, positions[0], positions[1])
					rep_node += (",%s" % total)
					coverage_e_s[element + "_" + strand] = float(total)
			to_add = classify(coverage_e_s, parameters[5], parameters[6])
			print(rep_node + "," + ",".join(to_add))
			rep_out.write(rep_node + "," + ",".join(to_add) + "\n")
	rep_out.close()
	print("\n(3) Done")
	
	print("\n(4) Sorting sequences")
	id_report = open(output_dir + "id_report.csv", "r").readlines()
	sorting_finished = sortSeqs(elements, output_dir, mapping_file, id_report)
	print("\n(4) Sorting sequences --- complete") if sorting_finished == "y" else print("\n(4) Sorting sequences --- incompleted")
	if mode != "-t" : os.system("rm -r " + output_dir + "/Temp")
	print("\n(4) Done")
		
if __name__ == "__main__":
	main()