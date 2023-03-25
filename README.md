# str-ngs-simulator

str-ngs-simulator is a tool to simulate Next Generation Sequencing datasets that include stutter errors characteristic of STRs. 

## Installation

With conda

```(insert conda installation)
```

## Usage

Steps to run tool
```(insert usage)
1. Download external packages "pyfaidx" and "art"
2. Run __main__.py file with parameters "--coords", "--ref", "--art"
3. Confirm that output "test_dir" was created with output files
```
ART Download:
https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm

Example run
```
General format:
python __main__.py --coords <path to coordinate file> --ref <path to reference genome> --art <path to art package>

With local parameters:
python __main__.py --coords ref_HTT.bed --ref /storage/resources/dbase/human/hg19/hg19.fa --art /storage/ashen/NGS_simulator/art_illumina --output_dir test_set_out
```

Output format
- The output directory contains a few different files and folders:
- "combined1.fq": combined final simulated reads in forward orientation
- "combined2.fq": combined final simulated reads in reverse orientation
- "fasta": directory containing the reference sequence files used as templates for simulated data
- "fastq": directory containing the individual simulated reads for each group of stutter errors 


Please reference directory "example_files" for example input and output files
1. ref_HTT.bed: example coordinates file
2. test_set_out: example output directory
3. Use hg19.fa reference genome 



## Development Notes
Currently str-ngs-simulator is...
BUILD TEST #1 (11/21)

## Contact Us
If you have any questions about usage or issues with the tool, please reach create a Github issue or ask a qustion in the Discussion section.

## Contributing

(insert contributing)

## License
[MIT](https://choosealicense.com/licenses/mit/)
