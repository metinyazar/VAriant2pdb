import pandas as pd

def insert_df(filename):
    main_file =  pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_VEP.xlsx')
    pdb_file = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\SIFTS\\'+filename+'_uniprot_pdb.xlsx')
    uniprot_file = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\SIFTS\\'+filename+'_geneid_uniprot.xlsx')
    df4 = pd.concat([main_file[['Uploaded_variation']],pdb_file, uniprot_file, main_file[['GENE_ID', 'Feature','ENSP','Gene_symbol','Protein_change']]],axis=1)
    df4.to_excel ('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_VEP_uniprotid_pdb_0.xlsx',index=False)
