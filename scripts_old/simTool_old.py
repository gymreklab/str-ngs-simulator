

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
#output_dir contains three items: fasta fastq combined.fq 


import sys
import os
import shutil
import argparse

#Argparse
parser = argparse.ArgumentParser(description='Performs ART simulation including the characteristic stutter error to STRs')
parser.add_argument('add_prob', metavar='u', type=int, help='u being probability of adding additional copy of repeat')
args = parser.parse_args()

#Creates foldering structure
current = os.getcwd()
new_path = current + "/" + sys.argv[7]

shutil.rmtree(new_path)
os.mkdir(new_path)
os.mkdir(new_path + "/fasta")
os.mkdir(new_path + "/fastq")

#Reading in files
input_fasta = sys.argv[5]
input_repeat_coords = sys.argv[6]
with open(input_fasta,'r') as file:
    input_fasta_contents = file.read()
    input_fasta_contents_split = input_fasta_contents.split()
    #input_fasta_contents_split: item 0 is key, item 1 is sequence
with open(input_repeat_coords,'r') as file:
    input_repeat_coords_contents = file.read()

#Locate repeat sequence
#print input_fasta_contents
#print input_repeat_coords_contents
coords_list = input_repeat_coords_contents.split()
start = int(coords_list[1])
end = int(coords_list[2])
repeat_seq = coords_list[4]
preflank = input_fasta_contents_split[1][0:start]
repeat = input_fasta_contents_split[1][start:end]
postflank = input_fasta_contents_split[1][end:]
#print preflank + "\n"
#print repeat + "\n"
#print postflank + "\n"

#Calculate HipSTR error model for values 
delta = [-3,-2,-1,0,1,2,3] #can be 0, >0, <0

values_dict = {}
u = float(sys.argv[1])
d = float(sys.argv[2])
rho = float(sys.argv[3])

for delt in delta:
    if(delt == 0):
	values_dict[delt] = 1-u-d
    elif(delt > 0):
	values_dict[delt] = (u*rho*(1-rho)**(delt-1))
    elif(delt < 0):
	values_dict[delt] = d*rho*(1-rho)**(-delt-1)



#Only keep percentages greater than p_thresh
p_thresh = float(sys.argv[4])
values_edit = {k:v for (k,v) in values_dict.items() if v > p_thresh}
#for k,v in values_edit.iteritems():
#	print k,v


#Create fasta files for all sequences with 
for k,v in values_edit.iteritems():
    stutter = k
    current_repeat = repeat
    new_fa = ""
    new_fa = new_fa + input_fasta_contents_split[0] + "\n" + preflank
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
    new_fa += postflank
    file_name = "output" + str(stutter) + ".fa"
    save_path = new_path + "/fasta"
    #print save_path
    completeName = os.path.join(save_path, file_name)
    #print completeName
    f = open(completeName, "w")
    f.write(new_fa)
    f.close()

#Run ART based on percentages
    value = int(v*1000)
    #print value
    #number of times to run ART
    #ART parameters: ref (fa file), read_length (100), coverage (100), insert (350), sd (50), output (fq file name)
    for val in range(value):
	ref = completeName
	read_length = 100
	coverage = 100
	insert = 350
	sd = 50
	name = "output" + str(stutter) + "_" + str(val)
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
master_fq = ""
for subdir, dirs, files in os.walk(new_path + "/fastq"):
    for dir in dirs:
	current_path = new_path + "/fastq/" + dir
	file_name = "/" + dir[:-4] + "1.fq"
	current_path += file_name
	with open(current_path) as file:
    	    content = file.read()
	master_fq += content + "\n"
f = open("/storage/ashen/NGS_simulator/test_dir/combined.fq", "w")
f.write(master_fq)
f.close()
   


