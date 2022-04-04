import pandas as pd



def filter_pdb(filename): #pdb olanları ve olmayanları ayırmak için
    main_file = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename +'_VEP_uniprotid_pdb_3.xlsx')
    with_pdb = main_file.dropna(subset=['Mapped_PDB'])
    with_pdb.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_VEP_uniprotid_pdb_4_withpdb.xlsx',index=False)
    no_pdb = main_file[main_file['Mapped_PDB'].isnull()]
    no_pdb.drop(['Amino acid', 'DSSP_8state','Relative ASA'], axis=1, inplace=True)
    no_pdb.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\' + filename + '_VEP_uniprotid_pdb_4_nopdb.xlsx',index=False)

