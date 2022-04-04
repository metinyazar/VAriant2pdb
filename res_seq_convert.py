import pandas as pd
import warnings

warnings.filterwarnings("ignore")




def res_seq_features(filename):
    df = pd.read_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_VEP_uniprotid_pdb_1.xlsx')
    df['wt_res'] = df['Protein_change'].str[0]
    df['mut_res'] = df['Protein_change'].str[-1]
    polar = ['R','K','E','D','Q','N']
    neutral = ['G','A','S','T','P','H','Y']
    hydrophobic = ['C','V','L','I','M','F','W']
    for i in polar:
        df.loc[(df.wt_res == i), 'wt_res_h'] = 'Polar'
        df.loc[(df.mut_res == i), 'mut_res_h'] = 'Polar'
    for j in neutral:
        df.loc[(df.wt_res == j), 'wt_res_h'] = 'Neutral'
        df.loc[(df.mut_res == j), 'mut_res_h'] = 'Neutral'
    for k in hydrophobic:
        df.loc[(df.wt_res == k), 'wt_res_h'] = 'Hydrophobic'
        df.loc[(df.mut_res == k), 'mut_res_h'] = 'Hydrophobic'
    small = ['G','A','S','C','T','P','D']
    medium = ['N','V','E','Q','I','L']
    large = ['M','H','K','F','R','Y','W']
    for i in small:
        df.loc[(df.wt_res == i), 'wt_res_v'] = 'Small'
        df.loc[(df.mut_res == i), 'mut_res_v'] = 'Small'
    for j in medium:
        df.loc[(df.wt_res == j), 'wt_res_v'] = 'Medium'
        df.loc[(df.mut_res == j), 'mut_res_v'] = 'Medium'
    for k in large:
        df.loc[(df.wt_res == k), 'wt_res_v'] = 'Large'
        df.loc[(df.mut_res == k), 'mut_res_v'] = 'Large'
    low_polarizability = ['G','A','S','D','T']
    medium_polarizability = ['C','P','N','V','E','Q','I','L']
    high_polarizability = ['K','M','H','F','R','Y','W']
    for i in low_polarizability:
        df.loc[(df.wt_res == i), 'wt_res_z'] = 'Small_polarizability'
        df.loc[(df.mut_res == i), 'mut_res_z'] = 'Small_polarizability'
    for j in medium_polarizability:
        df.loc[(df.wt_res == j), 'wt_res_z'] = 'Medium_polarizability'
        df.loc[(df.mut_res == j), 'mut_res_z'] = 'Medium_polarizability'
    for k in high_polarizability:
        df.loc[(df.wt_res == k), 'wt_res_z'] = 'Large_polarizability'
        df.loc[(df.mut_res == k), 'mut_res_z'] = 'Large_polarizability'
    low_polar = ['L','I','F','W','C','M','V','Y']
    neutral_polar = ['P','A','T','G','S']
    high_polar = ['H','Q','R','K','N','E','D']
    for i in low_polar:
        df.loc[(df.wt_res == i), 'wt_res_p'] = 'Low_polar'
        df.loc[(df.mut_res == i), 'mut_res_p'] = 'Low_polar'
    for j in neutral_polar:
        df.loc[(df.wt_res == j), 'wt_res_p'] = 'Neutral_polar'
        df.loc[(df.mut_res == j), 'mut_res_p'] = 'Neutral_polar'
    for k in high_polar:
        df.loc[(df.wt_res == k), 'wt_res_p'] = 'High_polar'
        df.loc[(df.mut_res == k), 'mut_res_p'] = 'High_polar'
    acidic_f = ['D','E']
    basic_f = ['H','K','R']
    polar_f = ['C','G','N','Q','S','T','Y']
    nonpolar_f = ['A','F','I','L','M','P','V','W']
    for i in acidic_f:
        df.loc[(df.wt_res == i), 'wt_res_f'] = 'Acidic'
        df.loc[(df.mut_res == i), 'mut_res_f'] = 'Acidic'
    for j in basic_f:
        df.loc[(df.wt_res == j), 'wt_res_f'] = 'Basic'
        df.loc[(df.mut_res == j), 'mut_res_f'] = 'Basic'
    for k in polar_f:
        df.loc[(df.wt_res == k), 'wt_res_f'] = 'Polar'
        df.loc[(df.mut_res == k), 'mut_res_f'] = 'Polar'
    for m in nonpolar_f:
        df.loc[(df.wt_res == m), 'wt_res_f'] = 'Nonpolar'
        df.loc[(df.mut_res == m), 'mut_res_f'] = 'Nonpolar'
    acidic_s = ['D','E']
    basic_s = ['H','K','R']
    aromatic = ['F','W','Y']
    amide = ['N','Q']
    small_hydroxyl = ['S','T']
    sulfur_con = ['C','M']
    aliphatic = ['A','G','P','I','L','V']
    for i in acidic_s:
        df.loc[(df.wt_res == i), 'wt_res_s'] = 'Acidic'
        df.loc[(df.mut_res == i), 'mut_res_s'] = 'Acidic'
    for j in basic_s:
        df.loc[(df.wt_res == j), 'wt_res_s'] = 'Basic'
        df.loc[(df.mut_res == j), 'mut_res_s'] = 'Basic'
    for  k in aromatic:
        df.loc[(df.wt_res == k), 'wt_res_s'] = 'Aromatic'
        df.loc[(df.mut_res == k), 'mut_res_s'] = 'Aromatic'
    for m in amide:
        df.loc[(df.wt_res == m), 'wt_res_s'] = 'Amide'
        df.loc[(df.mut_res == m), 'mut_res_s'] = 'Amide'
    for x in small_hydroxyl:
        df.loc[(df.wt_res == x), 'wt_res_s'] = 'Small_hydroxyl'
        df.loc[(df.mut_res == x), 'mut_res_s'] = 'Small_hydroxyl'
    for y in sulfur_con:
        df.loc[(df.wt_res == y), 'wt_res_s'] = 'Sulfur_containing'
        df.loc[(df.mut_res == y), 'mut_res_s'] = 'Sulfur_containing'
    for z in aliphatic:
        df.loc[(df.wt_res == z), 'wt_res_s'] = 'Aliphatic'
        df.loc[(df.mut_res == z), 'mut_res_s'] = 'Aliphatic'

    df.to_excel('D:\\tez_onerisi\\SNP\\work folder\\predicting_stability_structural_effect\\uniprot\\'+filename+'_VEP_uniprotid_pdb_2.xlsx',index=False)