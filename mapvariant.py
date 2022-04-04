import numpy as np
import pandas as pd
import warnings
import os
import re
from prody import fetchPDB
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.DSSP import DSSP
warnings.filterwarnings("ignore")



# ONE LETTER CODES
one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q', 'ASP': 'D',
              'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y', 'ARG': 'R',
              'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A', 'GLY': 'G',
              'PRO': 'P', 'CYS': 'C', ' DG': 'GTP', ' DA': 'ATP', ' DT': 'TTP', ' DC':'CTP', ' DI':'ITP',
              ' DU':'UTP', 'UNK': 'UNK'}

def is_het(residue):
    res = residue.id[0]
    if res != ' ':
        return True
    # if res == "W": # Just for water
    #     return True
    else:
        return False


def get_one_letter_aa_identifier(pdb_name):
    if len(pdb_name) == 4:
        pdbl = PDBList()
        pdb_id_list = []
        hetero_atoms = []
        parser = PDBParser(PERMISSIVE=1)  # PERMISSIVE was set to 2 (silence warnings)
        structure = parser.get_structure('pdb_name', fetchPDB(pdb_name, compressed=False))
        pdb_header = (structure.header['name'])
        for model in structure:
            for chain in model:
                for residue in chain:
                    if not is_het(residue):
                        ResId = residue.get_id()[1]
                        ResName = residue.get_resname()
                        pdb_id_list.append(one_letter[ResName] + str(ResId))
                    if is_het(residue):
                        hetero_atoms.append(residue.get_resname() + str(residue.get_id()[1]))

    return pdb_id_list

def getchain_id(pdb_id):
    pdbl = PDBList()
    parser = PDBParser(PERMISSIVE=1)  # PERMISSIVE was set to 2 (silence warnings)
    pdb = parser.get_structure(pdb_id, fetchPDB(pdb_id, compressed=False))
    data = []
    for chain in pdb.get_chains():
        data.append(chain.id)
    return data

def map_variant_pdb(filename):
    df = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_VEP_uniprotid_pdb_0.xlsx')
    index = [2]
    data2 = []
    data4 = []
    for k in index:
        pdbs = list(df['PDB'].iloc[k].split(';'))
        protein_change = list(map(lambda i: i[:-1], df['Protein_change'].iloc[k].split(',')))
        data = []
        data3 = []
        for j in pdbs:
            try:
                resid = get_one_letter_aa_identifier(j)
                a = list(set(protein_change) & set(resid))
                if len(a) == 0:
                    myfile = j + ".pdb"
                    ## If file exists, delete it ##
                    if os.path.isfile(myfile):
                        os.remove(myfile)
                        print("Deleted", j)
                    else:  ## Show an error ##
                        print("Error: %s file not found" % myfile)
                else:
                    data.append(j)
                    data3.append(a)
                    break #=>eğer sadece tek bir pdb dosyası yeterliyse break aktif edilir
            except Exception as Err:
                print("Error in", j)
        data2.append(data)
        data4.append(data3)
        print(k)
    df2 = pd.DataFrame()
    df3 = pd.DataFrame()
    df2['PDB'] = pd.Series(data2)
    df3['mapped_residue'] = pd.Series(data4)
    df5 = pd.concat([df2,df3],axis=1)
    df5.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_aradosya_pdb.xlsx',index=False)

map_variant_pdb('cancer_vusex')