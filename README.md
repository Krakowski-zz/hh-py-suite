## HH-py-suite
This project updates addss.pl and reformat.pl script from [hh-suite package](https://github.com/soedinglab/hh-suite).
The main goal is to make the scripts compatible with new versions of [psipred](https://github.com/psipred/psipred) (4.0) and 
[dssp](https://github.com/PDB-REDO/dssp) (4.0) and in result support mmCif format of the spatial protein structure.

The scripts were developed by Kamil Krakowski, Anna Semik and Natalia Rutecka as a part 
of "Architecture of large bioinformatics projects" course.

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
sequences in a way that will unable to later perform a sensitive sequence searching.

## Installation and requirements
The pipeline was written in Python 3.10 and uses the libraries specified in 
the requirements.txt file. 

Before you run the pipeline you also need to install:
1. [HHSuite](https://github.com/soedinglab/hh-suite)
2. [dssp](https://github.com/PDB-REDO/dssp)
3. [psipred](https://github.com/psipred/psipred)

## Parameters and configuration file

## Typical usage examples

## Examplary files
Examples of an input fasta file and an output file can be found in "Data" directory

