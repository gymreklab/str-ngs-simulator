#!/bin/bash

#input variables: u, d, rho, coverage

u=$1
d=$2
rho=$3
coverage=$4
ref="/storage/resources/dbase/human/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
coords="ref_HTT.bed"

#Runs simTool.py using the chosen parameters
python /storage/ashen/NGS_simulator/str-ngs-simulator/simTool.py --u $u --d $d --rho $rho --coverage $coverage --coords $coords --ref $ref >/dev/null

echo "INPUTS"
echo "u = " $u
echo "d = " $d
echo "rho = " $rho
echo "coverage = " $coverage
echo

total_reads=0
reads=0
expected_reads=0
prob=0
expected_prob=0
difference=0

#Calculate total reads
for dir in $(ls /storage/ashen/NGS_simulator/test_dir/fastq)
do
    dir_name=${dir%????}
    reads=$(wc -l < /storage/ashen/NGS_simulator/test_dir/fastq/${dir_name}_dir/${dir_name}1.fq)
    total_reads=$(($total_reads + $reads))
done

#Calculate rest of the values
#calculate the HipSTR error values for normalization of probability distributions
allProbs=()
for dir in $(ls /storage/ashen/NGS_simulator/test_dir/fastq)
do
    dir_name=${dir%????}
    dir_error="${dir_name:6}"
    reads=$(wc -l < /storage/ashen/NGS_simulator/test_dir/fastq/${dir_name}_dir/${dir_name}1.fq)
    prob=$(($reads / $total_reads))
    prob_dec=$(echo "scale=5; $reads / $total_reads" | bc)
#Application of HipSTR error model 
    if [ $dir_error == 0 ]
    then
        expected_prob=`echo "scale=5; 1 - $u - $d" | bc`
    elif [ $dir_error -gt 0 ]
    then
        inverse_rho=`echo "1 - $rho" | bc`
        error_1=`echo "$dir_error - 1" | bc`
        exp_term=`echo "scale=5; e($error_1*l($inverse_rho))" | bc -l`
        expected_prob=`echo "scale=5; $u * $rho * $exp_term" | bc`
    elif [ $dir_error -lt 0 ]
    then
        inverse_rho=`echo "1 - $rho" | bc`
        error_1=`echo "$dir_error * -1" | bc`
        error_2=`echo "$error_1 - 1" | bc`
        exp_term=`echo "scale=5; e($error_2*l($inverse_rho))" | bc -l`
        expected_prob=`echo "scale=5; $d * $rho * $exp_term" | bc`
    fi
    allProbs+=($expected_prob)
done
totalProb=0
for value in "${allProbs[@]}"
do
     totalProb=`echo "$totalProb + $value" | bc`
done
echo "TotalProb" is $totalProb

for dir in $(ls /storage/ashen/NGS_simulator/test_dir/fastq)
do
    dir_name=${dir%????}
    dir_error="${dir_name:6}"
    reads=$(wc -l < /storage/ashen/NGS_simulator/test_dir/fastq/${dir_name}_dir/${dir_name}1.fq)
    prob=$(($reads / $total_reads))i
#echo $var | awk '{print int($1+0.5)}' 
    prob_dec=$(echo "scale=15; $reads / $total_reads" | bc)
#Application of HipSTR error model 
    if [ $dir_error == 0 ]
    then
	expected_prob=`echo "scale=5; 1 - $u - $d" | bc`
    elif [ $dir_error -gt 0 ]
    then
	inverse_rho=`echo "1 - $rho" | bc`
	error_1=`echo "$dir_error - 1" | bc`
	exp_term=`echo "scale=5; e($error_1*l($inverse_rho))" | bc -l`
	expected_prob=`echo "scale=5; $u * $rho * $exp_term" | bc`
    elif [ $dir_error -lt 0 ]
    then
        inverse_rho=`echo "1 - $rho" | bc`
	error_1=`echo "$dir_error * -1" | bc`
	error_2=`echo "$error_1 - 1" | bc`
	exp_term=`echo "scale=5; e($error_2*l($inverse_rho))" | bc -l`
        expected_prob=`echo "scale=5; $d * $rho * $exp_term" | bc`
    fi
    expected_prob=`echo "scale=15; $expected_prob / $totalProb" | bc`
    expected_reads_dec=`echo "$total_reads * $expected_prob" | bc`

    expected_reads=${expected_reads_dec%.*}    

    difference_1=`echo "$reads - $expected_reads" | bc`
    difference=`echo "scale=5; $difference_1 / $expected_reads" | bc`

    echo $dir_name
    echo "Actual Reads = " $reads
    echo "Expected Reads = " $expected_reads " " $expected_reads_dec
    echo "Actual Probability = " $prob_dec
    echo "Expected Probability = " $expected_prob
    echo "Difference = " $difference
    echo
done

#for value in "${allProbs[@]}"
#do
#     echo $value
#done

