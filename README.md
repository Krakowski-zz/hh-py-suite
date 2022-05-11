# hh-py-suite
This project aims at updating addss.pl script from hh-suite package (https://github.com/soedinglab/hh-suite) in few steps:
<li>rewrite its functionality to python</li>
<li>change blastpgp (from legacy-blast package) to psiblast (from ncbi-blast+ package)</li>
<li>make it compatible with psipred v.4.0 (https://github.com/psipred/psipred)</li>
<li>update it to use dssp 4.0 with support to mmCIF (https://github.com/PDB-REDO/dssp)
<li>prepare parsers for different requested outputs and inputs</li>



# Before you run the pipeline, you need to install: 
1. HHSuite: `wget https://github.com/soedinglab/hh-suite/releases/download/v3.3.0/hhsuite-3.3.0-SSE2-Linux.tar.gz; tar xvfz hhsuite-3.3.0-SSE2-Linux.tar.gz; export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"`
2. ...
