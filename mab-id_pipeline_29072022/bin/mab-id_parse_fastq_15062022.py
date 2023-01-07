from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
from textwrap import dedent
import gzip 
import numpy
import re
import argparse
import os.path

parser = argparse.ArgumentParser(prog='AbID-UMI parser',formatter_class=argparse.RawDescriptionHelpFormatter,description=dedent('''\
	         This will find, clean and add UMI's to the output of cutadapt.
		 Every read will have a class assigned to it (in the order of discovery/confidence):
	         --------------------------------
	             1	both parts of the SBC are found at the expected place
	             2	both parts of the SBC are found, but with a 1-bp InDel in the first UMI
	             3	the first (a) or second (b) SBC part is found at the expected place
	             4	the first (a) or second (b) SBC part, but with a 1-bp InDel in the first UMI
	             5	the second SBC part is found somewhere in the matched adapter
	             6	the first SBC part is found somewhere in the matched adapter
	             7	the constant sequence is found somewhere in the matched adapter
	         --------------------------------
		 Classes 3 through 6 are used when:
		     - the other SBC-part is mutated 
		     - the middle UMI is not of the expected length.
		 Classes 7 is used when both SBC-parts have been mutated
		 Reads failing all these classes are removed.
	         '''))

parser.add_argument("input_file", help="An output from cutadapt: i.e., ABBC_001.SBC_011.fq")
parser.add_argument("adapter_file", help="Adapter-fasta given to cutadapt")
parser.add_argument("output_file", help="The name of the output-file to write to.")

args = parser.parse_args()

# variables
input_file = args.input_file
adapter_file = args.adapter_file 
expected_adapter_len = 30
handle = gzip.open(args.output_file, "wt")

# get SBC-adapter-sequences
adap_list = []
this_name = re.sub(".fq", "", os.path.basename(input_file)).split('.')
this_name = ".".join(this_name[1:])
for fa_t, fa_s in SimpleFastaParser(open(adapter_file, "rt")) :
	if fa_t != this_name :
		continue
	else :
		expected_adapter_len = len(re.sub("[\$\^\>]","", fa_s))
		ss = fa_s.split("NNN")
		ss = ss[1:3]
		ss.append("AGGGCCGC")
		ss.append(ss[1][len(ss[0])-len(ss[0]):len(ss[0])])
		bye = ss.pop(1)
		adap_list = [ss[i] for i in[0,2,1]]

adap_list = [each_string.lower() for each_string in adap_list]

# parse reads
for title, seq, qual in FastqGeneralIterator(open(input_file, "rt")) :
	if seq[3:7] == adap_list[0] and seq[10:14] == adap_list[1] :
		UMI = (seq[0:3] + seq[7:10]).upper()
		title = title.replace(" ", "_")+"_class1:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[2:6] == adap_list[0] and seq[9:13] == adap_list[1] :
		UMI = (seq[0:2] + seq[6:9]).upper()
		title = title.replace(" ", "_")+"_class2:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[4:8] == adap_list[0] and seq[11:15] == adap_list[1] :
		UMI = (seq[1:4] + seq[8:11]).upper()
		title = title.replace(" ", "_")+"_class2:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[10:14] == adap_list[1] :
		UMI = (seq[0:3] + seq[7:10]).upper()
		title = title.replace(" ", "_")+"_class3b:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[9:13] == adap_list[1] :
		UMI = (seq[0:2] + seq[6:9]).upper()
		UMI = re.sub('0','N',UMI.zfill(6,))
		title = title.replace(" ", "_")+"_class4b:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[11:15] == adap_list[1] :
		UMI = (seq[1:4] + seq[8:11]).upper()
		title = title.replace(" ", "_")+"_class4b:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[3:7] == adap_list[0] :
		UMI = (seq[0:3] + seq[7:10]).upper()
		title = title.replace(" ", "_")+"_class3a:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[2:6] == adap_list[0] :
		UMI = (seq[0:2] + seq[6:9]).upper()
		UMI = re.sub('0','N',UMI.zfill(6,))
		title = title.replace(" ", "_")+"_class4a:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif seq[4:8] == adap_list[0] :
		UMI = (seq[1:4] + seq[8:11]).upper()
		title = title.replace(" ", "_")+"_class4a:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif bool(re.search(adap_list[1],seq)) :
		sbc2_start = seq.find(adap_list[1])
		UMI = (seq[numpy.clip(sbc2_start-10,0,expected_adapter_len):numpy.clip(sbc2_start-7,0,expected_adapter_len)] + seq[sbc2_start-3:sbc2_start]).upper()
		if len(UMI) != 6 :
			 UMI = re.sub('0','N',UMI.zfill(6,))
		title = title.replace(" ", "_")+"_class5:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif bool(re.search(adap_list[0],seq)) :
		sbc1_start = seq.find(adap_list[0])
		UMI = (seq[sbc1_start-3:sbc1_start] + seq[sbc1_start+4:sbc1_start+7]).upper()
		if len(UMI) != 6 :
			 UMI = re.sub('0','N',UMI.zfill(6,))
		title = title.replace(" ", "_")+"_class6:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	elif bool(re.search(adap_list[2],seq)) :
		const_start = seq.find(adap_list[2])
		UMI = seq[:const_start-4].upper()[:3] + seq[:const_start-4].upper()[-3:]
		if len(UMI) != 6 :
			 UMI = re.sub('0','N',UMI.zfill(6,))
		title = title.replace(" ", "_")+"_class7:UMI:" + UMI
		qual = qual[min([i for i, a in enumerate(seq) if a.isupper()]):]
		seq = seq[min([i for i, a in enumerate(seq) if a.isupper()]):]
		handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
	else :
		continue

handle.close()






