import pandas as pd


def uniprot_pdb(filename):
    ''' # Bu script her bir varyantın sifts_geneid_uniprot.py scriptinin çıktısı olan Uniprot ID'lerinin
           SIFTS => http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz dosyasındaki PDB'lerini
           bulmak içindir'''
    sifts_pdb = pd.read_csv('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\SIFTS\\uniprot_pdb.csv',header=1)
    vep_uniprot = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\SIFTS\\'+filename+'_geneid_uniprot.xlsx')
    vep_uniprot['UniProt_ID'] = vep_uniprot['UniProt_ID'].str.replace('\W', '')
    data2 = []
    #index2 = vep_uniprot.index.values
    index2= [0,1,2]
    for j in index2:
        if sifts_pdb['SP_PRIMARY'].str.contains(vep_uniprot['UniProt_ID'].iloc[j]).any() == True:
            b=sifts_pdb.index[sifts_pdb['SP_PRIMARY'].str.contains(vep_uniprot['UniProt_ID'].iloc[j])]
            data2.append(sifts_pdb['PDB'][b[0]])
        elif vep_uniprot['UniProt_ID'].iloc[j]== 'Notfound':
            data2.append('Not found')
        else:
            data2.append('Not found')
        print(j)
    df2 = pd.DataFrame(data2,columns=['PDB'])
    df2.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\SIFTS\\'+filename+'_uniprot_pdb.xlsx',index=False)
    return df2