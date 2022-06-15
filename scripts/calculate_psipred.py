"""
Functions used to calculate psipred secondary structure from the input fasta file
"""

import subprocess
import sys
import time
import shutil
import os


def parse_fasta(filein):
    """ Parser for .fasta file"""
    fin = open(filein, 'r')
    previous_line = fin.__next__()
    sequences = []
    name = ""
    sequence = ""
    for line in fin:
        if previous_line.strip()[0] == '>':
            if name != "":
                sequences.append((name, sequence))
            name = previous_line.strip()[1:]
            sequence = ""
        else:
            sequence += previous_line.strip()
        previous_line = line
    else:
        sequence += previous_line.strip()
        sequences.append((name, sequence))
    return sequences


def print_error_cont(err):
    print(''.join(['# ' + x for x in err.strip().split('\n')]))


def copy_input(input_path, input_copy):
    """ Copies the input file to a temporary path"""
    shutil.copy2(input_path, input_copy)


def reformat_input(input_path, reformat_script, input_copy, file_a3m, M):
    """ Reformats input to a3m format"""
    available_formats = ["fasta", "fas", "aln"]
    filename, ext = os.path.splitext(input_path)
    ext = ext[1:]

    if ext in available_formats:
        print(
            subprocess.check_output([reformat_script, ext, "a3m", input_copy, file_a3m, "-M ", M]).decode('ascii'))

    else:
        raise Exception(
            f"Invalid input extension. Accepted extensions are {available_formats}. Please add valid extension to your file.")


def parse_horiz(horiz_str):
    """ Parser for .horiz file"""
    sequence = ""
    predicted_structure = ""
    confidence = ""
    for part in horiz_str.split('\n\n'):
        part = part.strip()
        if not part or part.startswith("#"):
            continue
        components = part.split('\n')
        conf, pred, seq = components[0], components[1], components[2]
        sequence += seq[6:]
        predicted_structure += pred[6:]
        confidence += conf[6:]
    return sequence, '\n'.join([predicted_structure[x:x+100] for x in range(0, len(predicted_structure), 100)]), '\n'.join([confidence[x:x+100] for x in range(0, len(confidence), 100)])


def strip_alignment(filein, fileout):
    with open(fileout, 'w') as fout:
        for name, seq in parse_fasta(filein):
            fout.write('>' + name + '\n' + seq.translate({45: None}) + '\n')


def count_consensus(hhconsensus, file_a3m, file_consensus, file_a3m_cons):
    """ Builds consensus sequence using hhconsensus """
    print('- ' + str(time.strftime('%X %x')) + " building consensus sequence")
    cmd_hhconsensus = subprocess.Popen(
        [hhconsensus, "-i", file_a3m, "-s", file_consensus, "-o", file_a3m_cons],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_hhconsensus.communicate()
    if err:
        print('- ' + str(time.strftime('%X %x')) + " following problems occurred during consensus sequence building: ")
        print_error_cont(err.decode("utf-8"))


def make_blast_db(makeblastdb, file_fastaseqs, path_db):
    """ Creates blast database"""
    print('- ' + str(time.strftime('%X %x')) + " creating BLAST database")
    cmd_mkblastdb = subprocess.Popen(
        [makeblastdb, "-in", file_fastaseqs, "-dbtype", "prot", "-out", path_db],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_mkblastdb.communicate()
    if err:
        print('- ' + str(time.strftime('%X %x')) + "following problems occurred during BLAST database creation: ")
        print_error_cont(err.decode("utf-8"))


def count_pssm(psiblast, file_db, file_chk, file_consensus, file_psiblast):
    """ Runs psiblast to get a PSSM matrix"""
    print('- ' + str(time.strftime('%X %x')) + " running psiblast")
    cmd_psiblast = subprocess.Popen(
        [psiblast, "-db", file_db, "-out_pssm", file_chk, "-query", file_consensus,
         "-inclusion_ethresh", "0.001",
         "-num_descriptions", "5000", "-num_iterations", "3", "-num_alignments", "0", "-out", file_psiblast],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_psiblast.communicate()
    if err:
        print('- ' + str(time.strftime('%X %x')) + "following problems occurred when running psiblast: ")
        print_error_cont(err.decode("utf-8"))


def parse_pssm(chkparse, file_chk, file_mtx):
    """Parses PSSM matrix and writes it to a file"""
    print('- ' + str(time.strftime('%X %x')) + " parsing PSSM matrix")
    cmd_psiblast = subprocess.Popen([chkparse, file_chk], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
    output, err = cmd_psiblast.communicate()
    if not output:
        print('- ' + str(time.strftime('%X %x')) + "An error occurred when parsing PSSM")
        print(err.decode("utf-8"))
        sys.exit()
    open(file_mtx, "wb").write(output)


def run_psipred(psipred, w1, w2, w3, file_mtx, file_ss):
    """ Runs psipred"""
    print('- ' + str(time.strftime('%X %x')) + " running PSIPRED")
    cmd_psipred = subprocess.Popen(
        [psipred, file_mtx, w1, w2, w3],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_psipred.communicate()
    if not output:
        print('- ' + str(time.strftime('%X %x')) + "An error occurred when running PSIPRED")
        print(err.decode("utf-8"))
        sys.exit()
    open(file_ss, "wb").write(output)


def generate_horiz(psipass2, w2, file_psipred, file_ss):
    """Generates .horiz file"""
    print('- ' + str(time.strftime('%X %x')) + " generating .horiz file")
    cmd_psipass2 = subprocess.Popen(
        [psipass2, w2, "1", "1.0", "1.0", file_psipred, file_ss],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_psipass2.communicate()
    if not output:
        print('- ' + str(time.strftime('%X %x')) + "An error occurred when generating .horiz file")
        print(err.decode("utf-8"))
        sys.exit()
    result = output.decode("ASCII")
    sequence, predicted_structure, confidence = parse_horiz(result)
    return sequence, predicted_structure, confidence
