import pandas as pd


df_component = pd.read_csv('dataset/merge/df_combined_with_components.csv')
df_bsh = pd.read_csv('dataset/merge/BSH_compound_protein_bile_acid_detailed 12.csv')

protein_cols = [col for col in df_bsh.columns if col.startswith('Protein_')]
inhibitor_cols = [col for col in df_bsh.columns if col.startswith('Inhibitor_')]

df_protein = df_bsh[['BSH_Enzyme'] + protein_cols].drop_duplicates()
df_inhibitor = df_bsh[['Inhibitor'] + inhibitor_cols].drop_duplicates()
df_inhibitor = df_inhibitor[df_inhibitor['Inhibitor'] != 'VI']

properties = [
    'complexity',
    'h_bond_acceptor',
    'h_bond_donor',
    'rotatable_bonds',
    'molecular_weight',
    'polar_surface_area'
]

Bile_Acids = ['GCA', 'GCDCA', 'GDCA', 'GLCA', 'GUDCA', 'TCA', 'TCDCA', 'TDCA',
              'TLCA', 'TUDCA']

Bile_Acid_rows = []

for acid in Bile_Acids:
    row = {'Bile_Acid': acid}
    for prop in properties:
        col_name = f"{acid}_{prop}"
        row[prop] = df_bsh[col_name].iloc[0]  # Assuming single-row source
    Bile_Acid_rows.append(row)

df_bile_acid = pd.DataFrame(Bile_Acid_rows)


df_protein_inhibitor_bile_acid = df_component.copy()
df_protein_inhibitor_bile_acid = df_protein_inhibitor_bile_acid.merge(df_protein, on='BSH_Enzyme', how='left')
df_protein_inhibitor_bile_acid = df_protein_inhibitor_bile_acid.merge(df_inhibitor, on='Inhibitor', how='left')
df_protein_inhibitor_bile_acid = df_protein_inhibitor_bile_acid.merge(df_bile_acid, on='Bile_Acid', how='left')

df_protein_inhibitor_bile_acid.to_csv('dataset/generated/df_protein_inhibitor_bile_acid.csv')