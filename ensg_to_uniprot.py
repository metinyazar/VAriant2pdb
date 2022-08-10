import pandas as pd



def ensg_uniprot(ensg_geneid, ensemble_sift):
    ''' # This function is to find the UniProt IDs of each variant from the Gene IDs (ENSG encoded) from the
    Ensembl VEP system in the SIFTS file. (Can be downloaded from http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_ensembl.csv.gz.)'''
    vep_output = pd.read_excel(ensg_geneid)
    sifts_ensemble = pd.read_csv(ensemble_sift,header=1)
    data = []
    index = vep_output.index.values
    for i in index:
        if sifts_ensemble['GENE_ID'].str.contains(vep_output['GENE_ID'].iloc[i]).any() == True:
            a = sifts_ensemble.index[sifts_ensemble['GENE_ID'].str.contains(vep_output['GENE_ID'].iloc[i])]
            data.append(str([sifts_ensemble['SP_PRIMARY'].iloc[a[0]]]))
            print(i)
        else:
            data.append('Not found')
        print(i)
    df = pd.DataFrame(data, columns=['UniProt_ID'])
    df.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\SIFTS\\'+filename+'_geneid_uniprot.xlsx',index=False)
    return df


