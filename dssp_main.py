import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBList
import re
import warnings
from prody import fetchPDB

warnings.filterwarnings("ignore")

def dssp(filename):
    def dssp_res(pdb_code,res_number,chain_id):
        # initialize PDB downloader
        pdb_dl = PDBList()
        # download PDB files, save them in current directory
        pdb_dl.retrieve_pdb_file(pdb_code, pdir='./', file_format='pdb', overwrite=True)
        #parse structure
        p = PDBParser(PERMISSIVE=1)
        structure = p.get_structure(pdb_code, './pdb%s.ent' % pdb_code)
        # use only the first model
        model = structure[0]
        # calculate DSSP
        dssp = DSSP(model, './pdb%s.ent' % pdb_code,file_type='PDB',dssp='/usr/bin/dssp') # your dssp file location
        # extract sequence and secondary structure from the DSSP tuple
        sequence = ''
        sec_structure = ''
        for z in range(len(dssp)):
            a_key = list(dssp.keys())[z]
            sequence += dssp[a_key][1]
            sec_structure += dssp[a_key][2]
        # print extracted sequence and structure
        return list(dssp[chain_id,(' ', res_number, ' ')])
        #print(sequence)
        #print(sec_structure)
    def getchain_id(pdb_id):
        pdbl = PDBList()
        parser = PDBParser(PERMISSIVE=1)  # PERMISSIVE was set to 2 (silence warnings)
        pdb = parser.get_structure(pdb_id, fetchPDB(pdb_id, compressed=False))
        data = []
        for chain in pdb.get_chains():
            data.append(chain.id)
        return data

    df = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename +'_VEP_uniprotid_pdb_2.xlsx')
    index =df.index.values
    data2 = []
    for i in index:
        pdb_id = str(df['mapped_single_PDB'][i])
        data = []
        if pdb_id != 'nan':
            res_id1 = re.findall('\d+',str(df['Protein_change'][i]))
            res_id2 = [int(j) for j in res_id1]
            chain_id = getchain_id(pdb_id)
            for k in res_id2:
                for c in chain_id:
                    try:
                        if dssp_res(pdb_id,k, c)[1] == str(df['wt_res'][i]):
                            data.append(dssp_res(pdb_id,k, c))
                            break
                        else:
                            print("Error: dssp is not found residue %s " % k)
                    except KeyError as Err:
                        print('Error in',k)
                    except Exception as Err2:
                        data.append('PDB is not suitable for DSSP analysis')
                        break
        else:
            data.append('PDB is Not found')
        print(i)
        data2.append(data)
        df2 = pd.DataFrame()
        df2['DSSP'] = pd.Series(data2)
        df2.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_DSSP.xlsx',index=False)
