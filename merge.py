import pandas as pd


df_component = pd.read_csv('dataset/merge/df_combined_with_components.csv')
df_bsh = pd.read_csv('dataset/merge/BSH_compound_protein_bile_acid_detailed 12.csv')

protein_inhibitor_cols = [col for col in df_bsh.columns if col.startswith('Protein_') or col.startswith('Inhibitor_')]

df_protein_inhibitor= df_bsh[['BSH_Enzyme'] + protein_inhibitor_cols]

df_protein_inhibitor_merged = pd.merge(df_component, df_protein_inhibitor, on='BSH_Enzyme', how='left').drop_duplicates()

properties = [
    'complexity',
    'h_bond_acceptor',
    'h_bond_donor',
    'rotatable_bonds',
    'molecular_weight',
    'polar_surface_area'
]

prop_cols = [c for c in df_bsh.columns if any(c.endswith(f"_{prop}") for prop in properties)]
df_bile_acid_props = df_bsh[prop_cols].drop_duplicates()

df_protein_inhibitor_bile_acid = df_protein_inhibitor_merged.copy()

for prop in properties:
    def get_prop_value(row):
        bile_acid = row['Bile_Acid']
        col_name = f"{bile_acid}_{prop}"
        return df_bile_acid_props.at[0, col_name]

    df_protein_inhibitor_bile_acid[prop] = df_protein_inhibitor_merged.apply(get_prop_value, axis=1)

print('break')