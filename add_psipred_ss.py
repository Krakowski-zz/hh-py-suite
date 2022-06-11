"""
Todo: add description here

Authors: Kamil Krakowski, Natalia Rutecka, Anna Semik
"""
import argparse
import os
import shutil
import subprocess
import sys
import time


def read_conf_file(pth=os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__))) + "/config.cfg"):
    """ Reads a configuration file and returns a map from arguments to values"""
    config = {}
    with open(pth, 'r') as config_file:
        for line in config_file:
            nm, val = line.split('::')
            config[nm.strip()] = val.strip()
    return config


def parse_arguments():
    """ Parses arguments from config file and from command line.
    If there is a conflict, the command line parameter will be used."""
    config_args = read_conf_file()
    parser = argparse.ArgumentParser(description="")  # todo: add description here
    parser.add_argument("--input", required=True, help="Path to the input file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    requirements = ["hhconsensus", "reformat", "makeblastdb", "hhfilter", "psiblast", "chkparse", "psipred",
                    "psipred_weights", "psipred_weights2", "psipred_weights3", "psipred_weights_p2", "psipass2"]
    for req in requirements:
        parser.add_argument(f"--{req}", required=False, help=f"Path to the requirement {req}")

    parser.set_defaults(**config_args)
    for action in parser._actions:
        if action.dest in config_args:
            action.required = False
    parsed = parser.parse_args()
    return parsed


def set_temp_paths(input_path):
    """Returns dictionary containing temporary paths"""
    filename, ext = os.path.splitext(input_path)
    ext = ext[1:]
    dirname = "tmp/" + str(time.time_ns()) + "/"
    os.makedirs(os.path.dirname(dirname), exist_ok=True)
    filepath = dirname + filename
    tmp = {}
    tmp["input_copy"] = filepath + "." + ext
    tmp["file_fastaseqs"] = filepath + ".fastaseqs"
    tmp["file_a3m"] = filepath + ".a3m"
    tmp["file_a3m_with_consensus"] = filepath + "_cons.a3m"
    tmp["file_consensus"] = filepath + ".cons"
    tmp["file_chk"] = filepath + ".chk"
    tmp["file_mtx"] = filepath + ".mtx"
    tmp["file_ss"] = filepath + ".ss"
    tmp["file_psiblast_out"] = filepath + ".o_psiblast"
    tmp["file_psipred_out"] = filepath + ".o_psipred"
    tmp["db"] = filepath
    tmp["dir"] = dirname
    return tmp


def copy_input(input_path, tmp_dict):
    """ Copies the input file to a temporary path"""
    shutil.copy2(input_path, tmp_dict["input_copy"])


def reformat_input(input_path, reformat_script, tmp_dict):
    """ Reformats input to a3m format"""
    available_formats = ["fasta", "fas", "aln"]
    filename, ext = os.path.splitext(input_path)
    ext = ext[1:]

    if ext in available_formats:
        print(
            subprocess.check_output([reformat_script, ext, "a3m", tmp_dict["input_copy"], tmp_dict["file_a3m"], "-M first"]).decode(
                'ascii'))

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
        conf, pred, seq, ind = part.split('\n')
        sequence += seq[6:]
        predicted_structure += pred[6:]
        confidence += conf[6:]
    return sequence, predicted_structure, confidence


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


def strip_alignment(filein, fileout):
    with open(fileout, 'w') as fout:
        for name, seq in parse_fasta(filein):
            fout.write('>' + name + '\n' + seq.translate({45: None}) + '\n')


def print_error_cont(err):
    print(''.join(['# ' + x for x in err.strip().split('\n')]))


def count_consensus(hhconsensus, tmp):
    """ Builds consensus sequence using hhconsensus """
    print('- ' + str(time.strftime('%X %x')) + " building consensus sequence")
    cmd_hhconsensus = subprocess.Popen(
        [hhconsensus, "-i", tmp['file_a3m'], "-s", tmp['file_consensus'], "-o",
         tmp['file_a3m_with_consensus']],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_hhconsensus.communicate()
    if err:
        print('- ' + str(time.strftime('%X %x')) + " following problems occurred during consensus sequence building: ")
        print_error_cont(err.decode("utf-8"))


def make_blast_db(makeblastdb, tmp):
    """ Creates blast database"""
    print('- ' + str(time.strftime('%X %x')) + " creating BLAST database")
    cmd_mkblastdb = subprocess.Popen(
        [makeblastdb, "-in", tmp["file_fastaseqs"], "-dbtype", "prot", "-out", tmp["db"]],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_mkblastdb.communicate()
    if err:
        print('- ' + str(time.strftime('%X %x')) + "following problems occurred during BLAST database creation: ")
        print_error_cont(err.decode("utf-8"))


def count_pssm(psiblast, tmp):
    """ Runs psiblast to get a PSSM matrix"""
    print('- ' + str(time.strftime('%X %x')) + " running psiblast")
    cmd_psiblast = subprocess.Popen(
        [psiblast, "-db", tmp["db"], "-out_pssm", tmp["file_chk"], "-query", tmp["file_consensus"],
         "-inclusion_ethresh", "0.001",
         "-num_descriptions", "5000", "-num_iterations", "3", "-num_alignments", "0", "-out", tmp["file_psiblast_out"]],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_psiblast.communicate()
    if err:
        print('- ' + str(time.strftime('%X %x')) + "following problems occurred when running psiblast: ")
        print_error_cont(err.decode("utf-8"))


def parse_pssm(chkparse, tmp):
    """Parses PSSM matrix and writes it to a file"""
    print('- ' + str(time.strftime('%X %x')) + " parsing PSSM matrix")
    cmd_psiblast = subprocess.Popen([chkparse, tmp["file_chk"]], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
    output, err = cmd_psiblast.communicate()
    if not output:
        print('- ' + str(time.strftime('%X %x')) + "An error occurred when parsing PSSM")
        print(err.decode("utf-8"))
        sys.exit()
    open(tmp["file_mtx"], "wb").write(output)


def run_psipred(psipred, w1, w2, w3, tmp):
    """ Runs psipred"""
    print('- ' + str(time.strftime('%X %x')) + " running PSIPRED")
    cmd_psipred = subprocess.Popen(
        [psipred, tmp["file_mtx"], w1, w2, w3],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_psipred.communicate()
    if not output:
        print('- ' + str(time.strftime('%X %x')) + "An error occurred when running PSIPRED")
        print(err.decode("utf-8"))
        sys.exit()
    open(tmp["file_ss"], "wb").write(output)


def generate_horiz(psipass2, w2, tmp):
    """Generates .horiz file"""
    print('- ' + str(time.strftime('%X %x')) + " generating .horiz file")
    cmd_psipass2 = subprocess.Popen(
        [psipass2, w2, "1", "1.0", "1.0", tmp["file_psipred_out"], tmp["file_ss"]],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = cmd_psipass2.communicate()
    if not output:
        print('- ' + str(time.strftime('%X %x')) + "An error occurred when generating .horiz file")
        print(err.decode("utf-8"))
        sys.exit()
    result = output.decode("ASCII")
    sequence, predicted_structure, confidence = parse_horiz(result)
    # TODO: make sequence line length compatible with output format
    return sequence, predicted_structure, confidence


def save_output(outpath, sequence, predicted_structure, confidence, tmp):
    with open(outpath, 'w') as out_f:
        out_f.write(">Consensus\n")
        out_f.write(sequence + '\n')
        out_f.write(">ss_pred\n")
        out_f.write(predicted_structure + '\n')
        out_f.write(">ss_conf\n")
        out_f.write(confidence + '\n')
        out_f.write(open(tmp["file_a3m"]).read())
        print("Psipred secondary structure was successfully added to the multialignment file")


def main():
    args = parse_arguments()
    tmp = set_temp_paths(args.input)
    copy_input(args.input, tmp)
    reformat_input(tmp["input_copy"], args.reformat, tmp)
    count_consensus(args.hhconsensus, tmp)
    strip_alignment(tmp["input_copy"], tmp["file_fastaseqs"])
    make_blast_db(args.makeblastdb, tmp)
    count_pssm(args.psiblast, tmp)
    parse_pssm(args.chkparse, tmp)
    run_psipred(args.psipred, args.psipred_weights, args.psipred_weights2, args.psipred_weights3, tmp)
    seq, predicted_str, conf = generate_horiz(args.psipass2, args.psipred_weights2, tmp)
    save_output(args.output, seq, predicted_str, conf, tmp)
    # os.remove(tmp["dir"])


if __name__ == "__main__":
    main()
