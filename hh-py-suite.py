"""
Todo: add description here

Authors: Kamil Krakowski, Natalia Rutecka, Anna Semik
"""
import argparse
from calculate_psipred import *
from calculate_dssp import run_dssp


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
    parser.add_argument("--format", required=False, default="mmCif", help="Format of 3d structure file (pdb or mmCif)")
    requirements = ["hhconsensus", "reformat", "makeblastdb", "hhfilter", "psiblast", "chkparse", "psipred",
                    "psipred_weights", "psipred_weights2", "psipred_weights3", "psipred_weights_p2", "psipass2", "mkdssp"]
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
    tmp["dssp_file"] = filepath + ".dssp"
    return tmp


def save_output(outpath, sequence, predicted_structure, confidence, file_a3m, dssp_file):
    with open(outpath, 'w') as out_f:
        out_f.write(">Consensus\n")
        out_f.write(sequence + '\n')
        out_f.write(">ss_pred\n")
        out_f.write(predicted_structure + '\n')
        out_f.write(">ss_conf\n")
        out_f.write(confidence + '\n')
        out_f.write(open(file_a3m).read())
        print("Psipred secondary structure was successfully added to the multialignment file")
        out_f.write(open(dssp_file).read())
        print("Dssp file was succesfully added")

def main():
    args = parse_arguments()
    tmp = set_temp_paths(args.input)
    copy_input(args.input, tmp["input_copy"])
    reformat_input(args.input, args.reformat, tmp["input_copy"], tmp["file_a3m"])
    count_consensus(args.hhconsensus, tmp["file_a3m"], tmp["file_consensus"], tmp["file_a3m_with_consensus"])
    strip_alignment(tmp["input_copy"], tmp["file_fastaseqs"])
    make_blast_db(args.makeblastdb, tmp["file_fastaseqs"], tmp["db"])
    count_pssm(args.psiblast, tmp["db"], tmp["file_chk"], tmp["file_consensus"], tmp["file_psiblast_out"])
    parse_pssm(args.chkparse, tmp["file_chk"], tmp["file_mtx"])
    run_psipred(args.psipred, args.psipred_weights, args.psipred_weights2, args.psipred_weights3, tmp["file_mtx"], tmp["file_ss"])
    seq, predicted_str, conf = generate_horiz(args.psipass2, args.psipred_weights2, tmp["file_psipred_out"], tmp["file_ss"])
    run_dssp(args.mkdssp, tmp["file_consensus"], args.format, tmp["dir"], tmp["dssp_file"])
    save_output(args.output, seq, predicted_str, conf, tmp["file_a3m"], tmp["dssp_file"])
    # os.remove(tmp["dir"])


if __name__ == "__main__":
    main()
