"""
Functions used to calculate DSSP file from consensus file.
DSSP file contains secondary structure, geometrical features and solvent exposure of proteins.
"""

from Bio import SeqIO
from Bio.PDB.PDBList import PDBList
import subprocess
from calculate_psipred import print_error_cont
import time


def get_pdb_code(consensus_file):
    """Returns pdb code of a sequence from consensus file"""
    record = SeqIO.read(consensus_file, format="fasta")
    pdb_code = record.id.split("_")[0]
    return pdb_code

def download_pdb_structure(pdb_code, format, dir):
    """Downloads a 3d structure from PDB using pdb code. Supports MMCIF format"""
    pdb_handle = PDBList()
    pdb_handle.retrieve_pdb_file(pdb_code, file_format=format, pdir=dir)


def run_dssp(mkdssp, consensus_file, format_3d, pdb_dir, output_file):
    """ Calculates secondary structure, geometrical features and solvent exposure of proteins"""
    pdb_code = get_pdb_code(consensus_file)
    download_pdb_structure(pdb_code, format_3d, pdb_dir)

    if format_3d == "mmCif":
        pdb_path = f"{pdb_dir}{pdb_code.lower()}.cif"
    elif format_3d == "pdb":
        pdb_path = f"{pdb_dir}pdb{pdb_code.lower()}.ent"
    dssp_cmd = subprocess.Popen(
        [mkdssp, "-i", pdb_path, "-o", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = dssp_cmd.communicate()
    if err:
        print('- ' + str(
            time.strftime('%X %x')) + " following problems occurred during consensus sequence building: ")
        print_error_cont(err.decode("utf-8"))

