import shutil, sys, os, time, subprocess

def read_conf_file(pth="config.cfg"):
    config = {}
    with open(pth, 'r') as config_file:
        for line in config_file:
            nm, val = line.split('::')
            config[nm.strip()] = val.strip()
    return config

def parse_cmd_args():
    args = {}
    for x in range(len(sys.argv)):
        if sys.argv[x].startswith('-'):
            args[sys.argv[x][1:]] = sys.argv[x+1]
        else:
            continue
    return args

def convert_alignment(reformat, input_file, output_file, input_format, output_format, parameters=[]):
    print(subprocess.check_output([reformat, input_format, output_format, input_file, output_file]+parameters).decode('ascii'))

def parse_horiz(horiz_str):
    sequence = ""
    predicted_structure = ""
    confidence = ""
    for part in result.split('\n\n'):
        part = part.strip()
        if not part or part.startswith("#"):
            continue
        conf, pred, seq, ind = part.split('\n')
        sequence += seq[6:]
        predicted_structure += pred[6:]
        confidence += conf[6:]
    return (sequence, predicted_structure, confidence)

def print_error_cont(err):
    print (''.join(['# '+x for x in err.strip().split('\n')]))

required = ["input", "output", "hhconsensus", "reformat", "makeblastdb", "hhfilter", "psiblast", "chkparse", "psipred",
            "psipred_weights", "psipred_weights2", "psipred_weights3", "psipred_weights_p2", "psipass2"]
available_formats = ["fasta", "fas", "aln"]
config = read_conf_file()
param = parse_cmd_args()
tmp = {}

if not all([True if (x in param or x in config) else False for x in required ]):
    print("Required parameters were not provided\nFor more information type -help")
    sys.exit()


filename, ext = os.path.splitext(param["input"])
ext = ext[1:]
dirname = "tmp/"+str(time.time_ns())+"/"
os.makedirs(os.path.dirname(dirname), exist_ok=True)
filepath = dirname+filename

tmp["input_copy"] = filepath + "." + ext
tmp["file_a3m"] = filepath + ".a3m"
tmp["file_a3m_with_consensus"] = filepath + "_cons.a3m"
tmp["file_consensus"] = filepath + ".cons"
tmp["file_chk"] = filepath + ".chk"
tmp["file_mtx"] = filepath + ".mtx"
tmp["file_ss"] = filepath + ".ss"
tmp["file_psiblast_out"] = filepath + ".o_psiblast"
tmp["file_psipred_out"] = filepath + ".o_psipred"

#create copy of input file
shutil.copy2(param["input"], tmp["input_copy"])

#convert alignment file
if ext in available_formats:
    convert_alignment(config["reformat"], tmp["input_copy"], tmp['file_a3m'], ext, "a3m", ["-M first"])

#count consensus structure
print('- '+str(time.strftime('%X %x'))+" building consensus sequence")
cmd_hhconsensus = subprocess.Popen([config["hhconsensus"], "-i", tmp['file_a3m'], "-s", tmp['file_consensus'], "-o", tmp['file_a3m_with_consensus']],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, err = cmd_hhconsensus.communicate()
if err:
    print('- '+str(time.strftime('%X %x'))+" following problems occurred during consensus sequence building: ")
    print_error_cont(err.decode("utf-8"))


#make blast database
print('- '+str(time.strftime('%X %x'))+" creating BLAST database")
cmd_mkblastdb = subprocess.Popen([config["makeblastdb"], "-in", tmp['file_a3m'], "-dbtype", "prot", "-out", filepath],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, err = cmd_mkblastdb.communicate()
if err:
    print('- '+str(time.strftime('%X %x'))+"following problems occurred during BLAST database creation: ")
    print_error_cont(err.decode("utf-8"))

#count PSSM matrix
print('- '+str(time.strftime('%X %x'))+" running psiblast")
cmd_psiblast = subprocess.Popen([config["psiblast"], "-db", filepath, "-out_pssm", tmp["file_chk"], "-query", tmp["file_consensus"], "-inclusion_ethresh", "0.001",
                                "-num_descriptions", "5000", "-num_iterations", "3", "-num_alignments", "0", "-out", tmp["file_psiblast_out"]],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, err = cmd_psiblast.communicate()
if err:
    print('- '+str(time.strftime('%X %x'))+"following problems occurred when running psiblast: ")
    print_error_cont(err.decode("utf-8"))

#parse PSSM for PSIPRED
print('- '+str(time.strftime('%X %x'))+" parsing PSSM matrix")
cmd_psiblast = subprocess.Popen([config["chkparse"], tmp["file_chk"]], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, err = cmd_psiblast.communicate()
if not output:
    print('- '+str(time.strftime('%X %x'))+"An error occurred when parsing PSSM")
    print(err.decode("utf-8"))
    sys.exit()
open(tmp["file_mtx"], "wb").write(output)

#run PSIPRED
print('- '+str(time.strftime('%X %x'))+" running PSIPRED")
cmd_psipred = subprocess.Popen([config["psipred"], tmp["file_mtx"], config["psipred_weights"], config["psipred_weights2"], config["psipred_weights3"]],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, err = cmd_psipred.communicate()
if not output:
    print('- '+str(time.strftime('%X %x'))+"An error occurred when running PSIPRED")
    print(err.decode("utf-8"))
    sys.exit()
open(tmp["file_ss"], "wb").write(output)

#generate .horiz file
print('- '+str(time.strftime('%X %x'))+" generating .horiz file")
cmd_psipass2 = subprocess.Popen([config["psipass2"], config["psipred_weights_p2"], "1", "1.0", "1.0", tmp["file_psipred_out"], tmp["file_ss"]], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, err = cmd_psipass2.communicate()
if not output:
    print('- '+str(time.strftime('%X %x'))+"An error occurred when generating .horiz file")
    print(err.decode("utf-8"))
    sys.exit()
result = output.decode("ASCII")
sequence, predicted_structure, confidence = parse_horiz(result)

with open(param["output"], 'w') as out_f:
    out_f.write(">Consensus\n")
    out_f.write(sequence+'\n')
    out_f.write(">ss_pred\n")
    out_f.write(predicted_structure+'\n')
    out_f.write(">Confidence\n")
    out_f.write(confidence+'\n')
    out_f.write(open(tmp["file_a3m"]).read())


print("Psipred secondary structure was successfully added to the multialignment file")

#remove tmp files
print(dirname)
#os.remove(dirname)
