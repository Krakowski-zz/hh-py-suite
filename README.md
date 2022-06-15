## HH-py-suite
This project updates addss.pl and reformat.pl script from [hh-suite package](https://github.com/soedinglab/hh-suite).
The main goal is to make the scripts compatible with new versions of [psipred](https://github.com/psipred/psipred) (4.0) and 
[dssp](https://github.com/PDB-REDO/dssp) (4.0) and in result support mmCif format of the spatial protein structure.

The scripts were developed by Kamil Krakowski, Anna Semik and Natalia Rutecka as a part 
of the "Architecture of large projects in bioinformatics" course.

## Introduction
HH-suite is a software package for highly sensitive sequence searching and sequence alignment.
It uses a concise representation of MSAs, called profile HMMS, which represent both
the query sequence and the database sequences. It can also use extra information about the
proteins, such as secondary structure or solvent accessibility. Because of that, hh-suite
often allows to make inferences from remotely homologous relationships, 
even when BLAST fail to do it.

Before the search methods can be used, the input fasta sequence or alignment file
needs to be preprocessed. This was initially done using scripts reformat.pl and
addss.pl from HH-Suite package. Unfortunately, the scripts haven't been updated for 3
years and are currently incompatible with most of the software that they use. Our goal 
is to restore their functionalities and create a pipeline that will allow to process
sequences in a way that will enable to later perform a sensitive sequence searching.

## Installation and requirements
The pipeline was written in Python 3.10 and uses the libraries specified in 
the requirements.txt file. 

Before you run the pipeline you also need to install:
1. [HHSuite](https://github.com/soedinglab/hh-suite)
2. [dssp](https://github.com/PDB-REDO/dssp)
3. [psipred](https://github.com/psipred/psipred)

You can also use conda environment:
conda install -c bioconda hhsuite
conda install -c salilab dssp
conda install -c biocore psipred

## Parameters and configuration file
The main program hh-py-suite.py has the following parameters: 

  --input INPUT         Path to the input file
  
  --output OUTPUT       Path to the output file
  
  --format FORMAT       Format of a 3d structure file (pdb or mmCif), default: mmCif
  
  --hhconsensus HHCONSENSUS
                        Path to the requirement hhconsensus
                        
  --reformat REFORMAT   Path to the requirement reformat 
  
  --makeblastdb MAKEBLASTDB
                        Path to the requirement makeblastdb 
                        
  --hhfilter HHFILTER   Path to the requirement hhfilter
  
  --psiblast PSIBLAST   Path to the requirement psiblast
  
  --chkparse CHKPARSE   Path to the requirement chkparse
  
  --psipred PSIPRED     Path to the requirement psipred
  
  --psipred_weights PSIPRED_WEIGHTS
                        Path to the requirement psipred_weights
                        
  --psipred_weights2 PSIPRED_WEIGHTS2
                        Path to the requirement psipred_weights2
                        
  --psipred_weights3 PSIPRED_WEIGHTS3
                        Path to the requirement psipred_weights3
                        
  --psipred_weights_p2 PSIPRED_WEIGHTS_P2
                        Path to the requirement psipred_weights_p2
                        
  --psipass2 PSIPASS2   Path to the requirement psipass2
  
  --mkdssp MKDSSP       Path to the requirement mkdssp
  
  The parameters can either be passed directly using command line or written into a config file. Each line of the file should contain name of the requirement and its path in a format: name:: path (see config file example in the data folder). If a parameter is specified both in a config file and via command line, the command line option will be used. 

## Typical usage
` python3 hh-py-suite.py --input sequences.fas --output sequences_extended.aln`
## Examplary files
Examples of an input fasta file, config file and output file can be found in the data directory.
