#!/usr/bin/env python

# Andrew Shen
# January 2021
# Performs ART simulation including the characteristic stutter error to STRs

# Main method for ease of packaging
def main():

    # Import packages
    import sys
    import os
    from os import path
    import shutil
    import argparse
    from pyfaidx import Fasta

    # Read in input parameters
    parser = argparse.ArgumentParser(description='Performs ART simulation including the characteristic stutter error to STRs')
    parser.add_argument("--u", help="Probability of adding additional copy of repeat", type=float, default=0.05)
    parser.add_argument("--d", help="Probability of deleting copy of repeat", type=float, default=0.05)
    parser.add_argument("--rho", help="Size of stutter-induced changes", type=float, default=0.9)
    parser.add_argument("--p_thresh", help="Cutoff percentage for meaningful values", type=float, default=0.01)
    parser.add_argument("--coverage", help="Size of coverage", type=int, default=1000)
    parser.add_argument("--read_length", help="Length of each read", type=int, default=100)
    parser.add_argument("--insert", help="(insert)", type=int, default=350)
    parser.add_argument("--sd", help="Value of standard deviation", type=int, default=50)
    parser.add_argument("--window", help="Size of window around sequence", type=int, default=1000)
    parser.add_argument("--coords", help="Path to file containing coordinates for desired sequence", type=str, required=True)
    parser.add_argument("--ref", help="Path to reference genome", type=str, required=True)
    parser.add_argument("--art", help="Path to ART simulator package", type=str, required=True)
    parser.add_argument("--output_dir", help="Name of output directory", type=str, default="test_dir")
    args = parser.parse_args()

    validArguments = True
    if(args.u < 0 or args.u > 1):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value between and 1 for --u.")
    if(args.d < 0 or args.d > 1):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value between and 1 for --d.")
    if(args.rho < 0 or args.rho > 1):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value between and 1 for --rho.")
    if(args.p_thresh < 0 or args.p_thresh > 1):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value between and 1 for --p_thresh.")
    if(args.coverage < 0):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value above 0 for --coverage.")
    if(args.read_length < 0):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value above 0 for --read_length.")
    if(args.insert < 0):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value above 0 for --insert.")
    if(args.sd < 0):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value above 0 for --sd.")
    if(args.window < 0):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid value above 0 for --window.")
    if(path.exists(args.coords) == False):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid path for --coords.")
    if(path.exists(args.ref) == False):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid path for --ref.")
    if(path.exists(args.art) == False):
        validArguments = False
        print("INVALID ARGUMENT: Please input a valid path for --art.")

    if(validArguments != True):
        exit()


    # Creates foldering structure
    current = os.getcwd()
    new_path = current + "/" + args.output_dir
    #If path to output directory exists, deletes directory, then creates a new directory
    if(os.path.isdir(new_path)): shutil.rmtree(new_path)
    #Creates the output directory and the fasta and fastq sub-directories
    os.mkdir(new_path)
    os.mkdir(new_path + "/fasta")
    os.mkdir(new_path + "/fastq")


    # Locate repeat sequence
    ref_genome = args.ref
    ref_coords = args.coords
    with open(ref_genome,'r') as file:
        ref_read = file.read()
    with open(ref_coords,'r') as file:
        ref_coords_read = file.read()
        ref_coords_split = ref_coords_read.split()

    if(ref_coords_split[0][:3] != "chr"):
        print('INVALID ARGUMENT: Please input a path to a --coords file in the format "chrom start end length repeat" with chrom in form "chr#".')
        exit()

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


    # More flexible version of HipSTR error model
    lowDelta = 0
    highDelta = 0

    values_dict = {}
    u = args.u
    d = args.d
    rho = args.rho
    p_thresh = args.p_thresh

    # Calculate lowDelta value
    inRange = True
    currDelta = 0
    while(inRange):
        calculatedValue = 0
        if(currDelta == 0):
            calculatedValue = 1-u-d
        elif(currDelta > 0):
            calculatedValue = (u*rho*(1-rho)**(currDelta-1))
        elif(currDelta < 0):
            calculatedValue = d*rho*(1-rho)**(-currDelta-1)    
        if(calculatedValue > p_thresh):
      	    currDelta = currDelta - 1
        else:
	        inRange = False
    lowDelta = currDelta

    # Calculate highDelta value
    inRange = True
    currDelta = 0
    while(inRange):
        calculatedValue = 0
        if(currDelta == 0):
            calculatedValue = 1-u-d
        elif(currDelta > 0):
            calculatedValue = (u*rho*(1-rho)**(currDelta-1))
        elif(currDelta < 0):
            calculatedValue = d*rho*(1-rho)**(-currDelta-1)
        if(calculatedValue > p_thresh):
            currDelta = currDelta + 1
        else:
            inRange = False
    highDelta = currDelta

    for delt in range(lowDelta+1, highDelta, 1):
        if(delt == 0):
           values_dict[delt] = 1-u-d
        elif(delt > 0):
           values_dict[delt] = (u*rho*(1-rho)**(delt-1))
        elif(delt < 0):
           values_dict[delt] = d*rho*(1-rho)**(-delt-1)

    # Only keep percentages greater than p_thresh
    values_edit = {k:v for (k,v) in values_dict.items() if v > p_thresh}

    # Create fasta files for all sequences
    art_path = args.art
    for k,v in values_edit.items():
        stutter = k
        header = ">" + str(chrom) + ":" + str(start) + "-" + str(end) + "_" + str(k)
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
            new_fa += current_repeat
            new_fa += str(postflank)
        file_name = "output" + str(stutter) + ".fa"
        save_path = new_path + "/fasta"
        completeName = os.path.join(save_path, file_name)
        f = open(completeName, "w")
        f.write(new_fa)
        f.close()

        # Run ART based on percentages
        value = int(v*int(args.coverage))
        ref = completeName
        read_length = args.read_length
        coverage = value
        insert = args.insert
        sd = args.sd
        name = "output" + str(stutter)
        dir_name = new_path + "/fastq/" + name
        os.mkdir(dir_name + "_dir")
        output = dir_name + "_dir/" + name
        cmd = art_path +' -sam -i '+ref+' -p -l '+str(read_length)+' -f '+str(coverage)+' -m '+str(insert)+' -s '+str(sd)+' -o '+output + " > /dev/null"
        os.system(cmd)

    # Combine all the fq files
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
    cwd = os.getcwd()
    f = open(new_path + "/combined1.fq", "w")
    f.write(master_fq_1)
    f.close()
    f2 = open(new_path+ "/combined2.fq", "w")
    f2.write(master_fq_2)
    f2.close()

    print("stutterSim ran successfully.")


if __name__ == "__main__":
    main()
