

#Andrew Shen
#January  2021

#Performs ART simulation including the characteristic stutter error to STRs
#To run: python simTool.py u d rho p_thresh input_fasta input_repeat_coords output_dir
	#sys.argv[1]: u being probability of adding additional copy of repeat (custom 0.05)
	#sys.argv[2]: d being probability of deleting copy of repeat (custom 0.05)
	#sys.argv[3]: rho being size of stutter-induced changes (custom 0.9)
	#sys.argv[4]: p_thresh being cutoff percentage 
	#sys.argv[5]: input_fasta being path to fasta file containing sequence
	#sys.argv[6]: input_repeat_coords being being path to file containing coordinates of tandem repeat
	#sys.argv[7]: output_dir being name of output directory
	#sys.argv[8]: ref_genome being path to the reference genome file
#output_dir contains three items: fasta fastq combined.fq 


import sys
import os
import shutil
import argparse
from pyfaidx import Fasta

#Argparse
parser = argparse.ArgumentParser(description='Performs ART simulation including the characteristic stutter error to STRs')
parser.add_argument("--u", help="Probability of adding additional copy of repeat", type=str, default=0.05)
parser.add_argument("--d", help="Probability of deleting copy of repeat", type=str, default=0.05)
parser.add_argument("--rho", help="Size of stutter-induced changes", type=str, default=0.9)
parser.add_argument("--p_thresh", help="Cutoff percentage for meaningful values", type=str, default=0.01)
parser.add_argument("--coverage", help="Size of coverage", type=str, default=1000)
parser.add_argument("--read_length", help="Length of each read", type=str, default=100)
parser.add_argument("--insert", help="(insert)", type=str, default=350)
parser.add_argument("--sd", help="Value of standard deviation", type=str, default=50)
parser.add_argument("--window", help="Size of window around sequence", type=str, default=1000)
parser.add_argument("--coords", help="Path to file containing coordinates for desired sequence", type=str, required=True)
parser.add_argument("--ref", help="Path to reference genome", type=str, required=True)
parser.add_argument("--output_dir", help="Name of output directory", type=str, default="test_dir")
args = parser.parse_args()

#Creates foldering structure
current = os.getcwd()
new_path = current + "/" + args.output_dir

shutil.rmtree(new_path)
os.mkdir(new_path)
os.mkdir(new_path + "/fasta")
os.mkdir(new_path + "/fastq")

#####
#Reading in files
#input_fasta = args.input_fasta
#input_repeat_coords = args.input_repeat_coords
#with open(input_fasta,'r') as file:
#    input_fasta_contents = file.read()
#    input_fasta_contents_split = input_fasta_contents.split()
#    #input_fasta_contents_split: item 0 is key, item 1 is sequence
#with open(input_repeat_coords,'r') as file:
#    input_repeat_coords_contents = file.read()

#####
#Locate repeat sequence
#print input_fasta_contents
#print input_repeat_coords_contents
#coords_list = input_repeat_coords_contents.split()
#start = int(coords_list[1])
#end = int(coords_list[2])
#repeat_seq = coords_list[4]
#preflank = input_fasta_contents_split[1][0:start]
#repeat = input_fasta_contents_split[1][start:end]
#postflank = input_fasta_contents_split[1][end:]
#print preflank + "\n"
#print repeat + "\n"
#print postflank + "\n"
ref_genome = args.ref
ref_coords = args.coords
with open(ref_genome,'r') as file:
    ref_read = file.read()
with open(ref_coords,'r') as file:
    ref_coords_read = file.read()
    ref_coords_split = ref_coords_read.split()

chrom = ref_coords_split[0]
start = int(ref_coords_split[1])
end = int(ref_coords_split[2])
repeat_length = int(ref_coords_split[3])
repeat_seq = ref_coords_split[4]
window = int(args.window)

ref_pyfaidx = Fasta(ref_genome)
region = ref_pyfaidx[chrom][start-window-1:end+window]
repeat = ref_pyfaidx[chrom][start-1:end]
preflank = ref_pyfaidx[chrom][start-window-1:start-1]
postflank = ref_pyfaidx[chrom][end:end+window]
#print repeat
#print 
#print preflank
#print
#print postflank

#Calculate HipSTR error model for values 
delta = [-3,-2,-1,0,1,2,3] #can be 0, >0, <0

values_dict = {}
u = float(args.u)
d = float(args.d)
rho = float(args.rho)

for delt in delta:
    if(delt == 0):
	values_dict[delt] = 1-u-d
    elif(delt > 0):
	values_dict[delt] = (u*rho*(1-rho)**(delt-1))
    elif(delt < 0):
	values_dict[delt] = d*rho*(1-rho)**(-delt-1)



#Only keep percentages greater than p_thresh
p_thresh = float(args.p_thresh)
values_edit = {k:v for (k,v) in values_dict.items() if v > p_thresh}
#for k,v in values_edit.iteritems():
#	print k,v

#print values_edit
#Create fasta files for all sequences with 
for k,v in values_edit.iteritems():
    stutter = k
    #FASTA header: >chr4:3073877-3075933
    header = ">" + str(chrom) + ":" + str(start) + "-" + str(end) + "_" + str(k)
    #print header
    current_repeat = str(repeat)
    new_fa = ""
    new_fa = new_fa + header + "\n" + str(preflank)
    if k == 0:
	new_fa += current_repeat
    elif k > 0:
	while k > 0:
	    current_repeat += repeat_seq
	    k-=1
	new_fa += current_repeat
    else:
	repeat_length = len(repeat_seq)
	subtract_size = k*repeat_length
	current_repeat = current_repeat[:subtract_size]
    new_fa += str(postflank)
    file_name = "output" + str(stutter) + ".fa"
    save_path = new_path + "/fasta"
    #print save_path
    completeName = os.path.join(save_path, file_name)
    #print completeName
    f = open(completeName, "w")
    f.write(new_fa)
    f.close()

#Run ART based on percentages
    #value = int(v*int(args.coverage))
    value = v*int(args.coverage)
    print value
    #number of times to run ART
    #ART parameters: ref (fa file), read_length (100), coverage (value based on percentages), insert (350), sd (50), output (fq file name)
    #for val in range(value):
    ref = completeName
    read_length = args.read_length
    coverage = value
    insert = args.insert
    sd = args.sd
    name = "output" + str(stutter)
    dir_name = new_path + "/fastq/" + name
    os.mkdir(dir_name + "_dir")
    output = dir_name + "_dir/" + name
    #print name
    #print output
    cmd='/storage/ashen/NGS_simulator/art_illumina -sam -i '+ref+' -p -l '+str(read_length)+' -f '+str(coverage)+' -m '+str(insert)+' -s '+str(sd)+' -o '+output
    os.system(cmd)
    #break
    #break

#Combine all the fq files
master_fq_1 = ""
master_fq_2 = ""
for subdir, dirs, files in os.walk(new_path + "/fastq"):
    for dir in dirs:
	current_path1 = new_path + "/fastq/" + dir
	current_path2 = new_path + "/fastq/" + dir
	file_name1 = "/" + dir[:-4] + "1.fq"
	file_name2 = "/" + dir[:-4] + "2.fq"
	current_path1 += file_name1
	current_path2 += file_name2
	with open(current_path1) as file:
    	    content1 = file.read()
	master_fq_1 += content1 + "\n"
	with open(current_path2) as file:
            content2 = file.read()
	master_fq_2 += content2 + "\n"
f = open("/storage/ashen/NGS_simulator/test_dir/combined1.fq", "w")
f.write(master_fq_1)
f.close()
f2 = open("/storage/ashen/NGS_simulator/test_dir/combined2.fq", "w")
f2.write(master_fq_2)
f2.close()

   


