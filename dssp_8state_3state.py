import pandas as pd


def state_conversion(filename):
    main_file = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename +'_VEP_uniprotid_pdb_4_withpdb.xlsx')
    main_file['DSSP_8state'] = main_file['DSSP_8state'].str.replace(r'\'','')
    main_file['DSSP_8state'] = main_file['DSSP_8state'].str.strip()
    coil = ['-','I', 'T', 'S']
    helix = ['G','H']
    extended = ['B','E']
    for i in coil:
        main_file.loc[(main_file['DSSP_8state'] == i),'DSSP_3state'] = 'C'
    for j in helix:
        main_file.loc[(main_file['DSSP_8state'] == j), 'DSSP_3state'] = 'H'
    for k in extended:
        main_file.loc[(main_file['DSSP_8state'] == k), 'DSSP_3state'] = 'E'
    main_file.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_VEP_uniprotid_pdb_5_withpdb.xlsx',index=False)
