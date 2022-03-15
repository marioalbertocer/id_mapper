import os, sys
from Bio import SeqIO

def setBlastInputs(mapping, elements, out_dir):
	tempDir = out_dir +  "Temp/"
	if not os.path.exists(tempDir): os.system("mkdir " + tempDir)
	concaFile = tempDir + "all_query.fasta"
	concaFile_w = open(concaFile, "w")
	assemblies = []
	for l in mapping:
		if "ref genome" in l.split(",")[0].strip():
			f, p, element, access = map(lambda x: x.strip(), l.split(",")[1:5])
			if element in elements:
				print("Query genome for element %s: %s" % (element,access))
				for rec in SeqIO.parse(p + "/" + f, "fasta"):
					concaFile_w.write(">%s_%s\n%s\n" % (element, access, rec.seq))
		if (l.split(",")[0]).strip() == "assembly":
			f, p, s = map(lambda x: x.strip(), l.split(",")[1:4])
			os.system("cp " + p + "/" + f + " " + tempDir)
			assemblies.append(tempDir + f + "," + s)
	return concaFile, assemblies

def runBlast(parameters, concaFile, assemblies):
	evalue, p_id, n_threads, out_dir, elements = parameters[:5]
	other_parameters = '-strand both -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"'
	blast_reports = []
	for assembly in assemblies:
		s = assembly.split(",")[1]
		assembly = assembly.split(",")[0]
		out_file = assembly.replace(".fasta", ".csv")
		os.system("makeblastdb -in " + assembly + " -dbtype nucl -parse_seqids")
		print('\nblastn -query %s -db %s -evalue %s -out %s -perc_identity %s -num_threads %s %s' % (concaFile, assembly, evalue, out_file, p_id, n_threads, other_parameters))
		os.system('blastn -query %s -db %s -evalue %s -out %s -perc_identity %s -num_threads %s %s' % (concaFile, assembly, evalue, out_file, p_id, n_threads, other_parameters))
		print(out_file + "-- Done")
		blast_reports.append("%s,%s" % (out_file, s))
	return blast_reports