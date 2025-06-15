import pandas as pd
import requests
import xml.etree.ElementTree as ET
from Bio.SeqUtils.ProtParam import ProteinAnalysis

'''
mass (g/L) = molarity (mol/L) × MW (g/mol)
1μM = (1×10)^−6 (mol/L)
mass = (1×10)^−6 (mol/L) × MW (g/mol)
mass = MW × 10^−6 (g/L)
1g = 10^6 (μg)
1L = 1000mL
1 (g/L) = 1000 (μg/mL)
mass = MW × 10^−6 (g/L) = MW × 10^−6 × 10^6 (μg/mL) = MW (μg/mL) --> means Concentration is 1 µM
μg/mL = μM × MW

100μM --> 100 * MW
200nM --> (200/1000) μM * MW = 0.2 * MW
'''
class MetadataGeneration:
    def __init__(self):
        self.bsh_data = None
        self.inhibitors_data = None
        self.bsh_enzymes_data = None
        self.bile_acid_data = None
        self.micromolar_enzyme_protein = 0.2
        self.micromolar_inhibitor = 100
        self.micromolar_bile_acids = 100

    def load_data(self):
        bsh_df = pd.read_csv('dataset/Clustering_data.csv')
        inhibitors_df = pd.read_csv('dataset/inhibitors.csv')
        bsh_enzymes_df = pd.read_csv('dataset/bsh_enzymes.csv')
        bile_acid_df = pd.read_csv('dataset/bile_acid.csv')
        self.inhibitors_data = inhibitors_df
        self.bsh_data = bsh_df
        self.bsh_enzymes_data = bsh_enzymes_df
        self.bile_acid_data = bile_acid_df

    def fetch_pubchem_compound(self, cid):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
        resp = requests.get(url)
        if not resp.ok:
            print(f"Failed for CID: {cid}")
            return None

        data = resp.json()
        compound = data.get("PC_Compounds", [{}])[0]

        props = {
            p['urn'].get('label', '') + ((" " + p['urn'].get('name', '')) if p['urn'].get('name') else ""): list(p['value'].values())[0]
            for p in compound.get('props', [])
        }

        core = {
            "cid": cid,
            "compound_canonicalized": props.get("Compound Canonicalized"),
            "complexity": props.get("Compound Complexity"),
            "h_bond_acceptor": props.get("Count Hydrogen Bond Acceptor"),
            "h_bond_donor": props.get("Count Hydrogen Bond Donor"),
            "rotatable_bonds": props.get("Count Rotatable Bond"),
            "fingerprint_substructure": props.get("Fingerprint SubStructure Keys"),
            "iupac_allowed": props.get("IUPAC Name Allowed"),
            "iupac_cas": props.get("IUPAC Name CAS-like Style"),
            "iupac_markup": props.get("IUPAC Name Markup"),
            "iupac_preferred": props.get("IUPAC Name Preferred"),
            "iupac_systematic": props.get("IUPAC Name Systematic"),
            "iupac_traditional": props.get("IUPAC Name Traditional"),
            "inchi": props.get("InChI Standard"),
            "inchikey": props.get("InChIKey Standard"),
            "logp": props.get("Log P XLogP3"),
            "mass_exact": props.get("Mass Exact"),
            "molecular_formula": props.get("Molecular Formula"),
            "molecular_weight": props.get("Molecular Weight"),
            "smiles_absolute": props.get("SMILES Absolute"),
            "smiles_canonical": props.get("SMILES Canonical"),
            "smiles_isomeric": props.get("SMILES Isomeric"),
            "polar_surface_area": props.get("Topological Polar Surface Area"),
            "weight_monoisotopic": props.get("Weight MonoIsotopic")
        }

        return core

    def get_data_details(self, source_df, id_column, fetch_function):
        results = []
        for item_id in source_df[id_column].dropna().unique():
            try:
                data = fetch_function(item_id)
                if data:
                    results.append(data)
            except Exception as e:
                print(f"Error fetching {item_id}: {e}")

        return pd.DataFrame(results)


    def get_enzyme_protein_details(self):
        result_dfs = {}
        accession_numbers = self.bsh_enzymes_data['PROTEIN_ACCESSION_NUMBER'].dropna().unique()

        for acc in accession_numbers:
            try:
                data = self.fetch_protein_xml_as_json(acc)
                df = pd.DataFrame([data])  # Wrap in list to form single-row DataFrame
                result_dfs[acc] = df
            except Exception as e:
                print(f"Error fetching {acc}: {e}")

        return result_dfs

    def fetch_protein_xml_as_json(self, accession):
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "protein",
            "id": accession,
            "retmode": "xml"
        }

        response = requests.get(url, params=params)
        root = ET.fromstring(response.text)

        # Extract sample fields
        docsum = {}
        for seq in root.findall(".//GBSeq"):
            docsum["accession"] = seq.findtext("GBSeq_primary-accession")
            docsum["definition"] = seq.findtext("GBSeq_definition")
            docsum["organism"] = seq.findtext("GBSeq_organism")
            docsum["sequence"] = seq.findtext("GBSeq_sequence")
        return docsum

    def analyze_protein(self, sequence):
        seq = sequence.upper().replace(" ", "").replace("\n", "")
        analysis = ProteinAnalysis(seq)
        return {
            "molecular_weight": round(analysis.molecular_weight(),3),
            "aromaticity": round(analysis.aromaticity(),3),
            "instability_index": round(analysis.instability_index(),3),
            "isoelectric_point": round(analysis.isoelectric_point(),3),
            "gravy": round(analysis.gravy(),3),
            "amino_acid_percent": {k: round(v, 3) for k, v in analysis.get_amino_acids_percent().items()}
        }

    def enrich_protein_df(self, protein_df):
        analysis_results = protein_df['sequence'].apply(self.analyze_protein)
        analysis_df = pd.json_normalize(analysis_results)
        return pd.concat([protein_df, analysis_df], axis=1)

    def generate_metadata(self):
        self.load_data()
        compounds = self.get_data_details(
            source_df=self.inhibitors_data,
            id_column='CID',
            fetch_function=self.fetch_pubchem_compound
        )
        protein_bsh = self.get_data_details(
            source_df=self.bsh_enzymes_data,
            id_column='PROTEIN_ACCESSION_NUMBER',
            fetch_function=self.fetch_protein_xml_as_json
        )
        bile_acids = self.get_data_details(
            source_df=self.bile_acid_data,
            id_column='CID',
            fetch_function=self.fetch_pubchem_compound
        )
        protein_bsh_enrich= self.enrich_protein_df(protein_bsh)
        compounds_detailed = pd.merge(compounds.add_prefix('inhibitor_'), self.inhibitors_data.add_prefix('inhibitor_'), how='left', left_on='inhibitor_cid', right_on='inhibitor_CID')
        BSH_compound_detailed = pd.merge(self.bsh_data.add_prefix('bsh_'), compounds_detailed,  how='left', left_on='bsh_Inhibitor', right_on='inhibitor_CODE')
        protein_bsh_detailed = pd.merge(protein_bsh_enrich.add_prefix('protein_'), self.bsh_enzymes_data, how='left', left_on='protein_accession', right_on='PROTEIN_ACCESSION_NUMBER')
        BSH_compound_protein_detailed = pd.merge(BSH_compound_detailed, protein_bsh_detailed, how='left', left_on='bsh_BSH enzyme', right_on='BSH_ENZYME')
        BSH_compound_protein_detailed['Inhibitor_Concentration'] = round(
                float(self.micromolar_inhibitor) *
                pd.to_numeric(BSH_compound_protein_detailed['inhibitor_molecular_weight'], errors='coerce')
        , 3)
        BSH_compound_protein_detailed['Enzyme_Protein_Concentration'] = round(
                float(self.micromolar_enzyme_protein) *
                pd.to_numeric(BSH_compound_protein_detailed['protein_molecular_weight'], errors='coerce')
        , 3)
        bile_acids_detailed = pd.merge(bile_acids, self.bile_acid_data, how='left', left_on='cid', right_on='CID')
        print(bile_acids_detailed['molecular_weight'])
        bile_acids_detailed['Bile_Acid_Concentration'] = (
                float(self.micromolar_bile_acids) *
                pd.to_numeric(bile_acids_detailed['molecular_weight'], errors='coerce')
        )
        BSH_compound_protein_detailed['TUDCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'TUDCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['TCDCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'TCDCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['TDCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'TDCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['TCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'TCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['TLCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'TLCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['GUDCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'GUDCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['GCDCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'GCDCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['GDCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'GDCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['GCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'GCA', 'Bile_Acid_Concentration'].values[0]
        BSH_compound_protein_detailed['GLCA_Concentration'] = bile_acids_detailed.loc[bile_acids_detailed['BILE_ACID'] == 'GLCA', 'Bile_Acid_Concentration'].values[0]


        original_columns = [
            'bsh_Run', 'bsh_BSH enzyme', 'Enzyme_Protein_Concentration', 'bsh_Inhibitor',
            'Inhibitor_Concentration', 'bsh_Replicate', 'bsh_Sample ID', 'TUDCA_Concentration',
            'TCA_Concentration', 'TLCA_Concentration', 'GUDCA_Concentration', 'GCDCA_Concentration',
            'GDCA_Concentration', 'GCA_Concentration', 'GLCA_Concentration', 'bsh_TUDCA', 'bsh_TCDCA',
            'bsh_TDCA', 'bsh_TCA', 'bsh_TLCA', 'bsh_GUDCA', 'bsh_GCDCA', 'bsh_GDCA', 'bsh_GCA',
            'bsh_GLCA', 'inhibitor_compound_canonicalized', 'inhibitor_complexity',
            'inhibitor_h_bond_acceptor', 'inhibitor_h_bond_donor', 'inhibitor_rotatable_bonds',
            'inhibitor_mass_exact', 'inhibitor_molecular_weight',
            'inhibitor_polar_surface_area', 'inhibitor_weight_monoisotopic', 'protein_organism',
            'protein_molecular_weight', 'protein_aromaticity', 'protein_instability_index',
            'protein_isoelectric_point', 'protein_gravy', 'protein_amino_acid_percent.A',
            'protein_amino_acid_percent.C', 'protein_amino_acid_percent.D', 'protein_amino_acid_percent.E',
            'protein_amino_acid_percent.F', 'protein_amino_acid_percent.G', 'protein_amino_acid_percent.H',
            'protein_amino_acid_percent.I', 'protein_amino_acid_percent.K', 'protein_amino_acid_percent.L',
            'protein_amino_acid_percent.M', 'protein_amino_acid_percent.N', 'protein_amino_acid_percent.P',
            'protein_amino_acid_percent.Q', 'protein_amino_acid_percent.R', 'protein_amino_acid_percent.S',
            'protein_amino_acid_percent.T', 'protein_amino_acid_percent.V', 'protein_amino_acid_percent.W',
            'protein_amino_acid_percent.Y'
        ]


        BSH_compound_protein_detailed = BSH_compound_protein_detailed[original_columns]

        BSH_compound_protein_detailed = BSH_compound_protein_detailed.rename(columns={
            'bsh_Run': 'Run',
            'bsh_BSH enzyme': 'BSH_Enzyme',
            'bsh_Inhibitor': 'Inhibitor',
            'bsh_Replicate': 'Replicate',
            'bsh_Sample ID': 'Sample_ID',
            'bsh_TUDCA': 'TUDCA',
            'bsh_TCDCA': 'TCDCA',
            'bsh_TDCA': 'TDCA',
            'bsh_TCA': 'TCA',
            'bsh_TLCA': 'TLCA',
            'bsh_GUDCA': 'GUDCA',
            'bsh_GCDCA': 'GCDCA',
            'bsh_GDCA': 'GDCA',
            'bsh_GCA': 'GCA',
            'bsh_GLCA': 'GLCA',
            'inhibitor_complexity': 'Inhibitor_Complexity',
            'inhibitor_h_bond_acceptor': 'Inhibitor_H_Bond_Acceptor',
            'inhibitor_h_bond_donor': 'Inhibitor_H_Bond_Donor',
            'inhibitor_rotatable_bonds': 'Inhibitor_Rotatable_Bonds',
            'inhibitor_mass_exact': 'Inhibitor_Mass_Exact',
            'inhibitor_molecular_weight': 'Inhibitor_Molecular_Weight',
            'inhibitor_polar_surface_area': 'Inhibitor_Polar_Surface_Area',
            'inhibitor_weight_monoisotopic': 'Inhibitor_Weight_Monoisotopic',
            'protein_organism': 'Protein_Organism',
            'protein_molecular_weight': 'Protein_Molecular_Weight',
            'protein_aromaticity': 'Protein_Aromaticity',
            'protein_instability_index': 'Protein_Instability_Index',
            'protein_isoelectric_point': 'Protein_Isoelectric_Point',
            'protein_gravy': 'Protein_GRAVY',
            'protein_amino_acid_percent.A': 'Protein_AA_Percent_A',
            'protein_amino_acid_percent.C': 'Protein_AA_Percent_C',
            'protein_amino_acid_percent.D': 'Protein_AA_Percent_D',
            'protein_amino_acid_percent.E': 'Protein_AA_Percent_E',
            'protein_amino_acid_percent.F': 'Protein_AA_Percent_F',
            'protein_amino_acid_percent.G': 'Protein_AA_Percent_G',
            'protein_amino_acid_percent.H': 'Protein_AA_Percent_H',
            'protein_amino_acid_percent.I': 'Protein_AA_Percent_I',
            'protein_amino_acid_percent.K': 'Protein_AA_Percent_K',
            'protein_amino_acid_percent.L': 'Protein_AA_Percent_L',
            'protein_amino_acid_percent.M': 'Protein_AA_Percent_M',
            'protein_amino_acid_percent.N': 'Protein_AA_Percent_N',
            'protein_amino_acid_percent.P': 'Protein_AA_Percent_P',
            'protein_amino_acid_percent.Q': 'Protein_AA_Percent_Q',
            'protein_amino_acid_percent.R': 'Protein_AA_Percent_R',
            'protein_amino_acid_percent.S': 'Protein_AA_Percent_S',
            'protein_amino_acid_percent.T': 'Protein_AA_Percent_T',
            'protein_amino_acid_percent.V': 'Protein_AA_Percent_V',
            'protein_amino_acid_percent.W': 'Protein_AA_Percent_W',
            'protein_amino_acid_percent.Y': 'Protein_AA_Percent_Y',
            'Inhibitor_Concentration': 'Inhibitor_IP_Conc',
            'Enzyme_Protein_Concentration': 'BSH_IP_Conc',
            'TUDCA_Concentration': 'TUDCA_Conc',
            'TCDCA_Concentration': 'TCDCA_Conc',
            'TDCA_Concentration': 'TDCA_Conc',
            'TCA_Concentration': 'TCA_Conc',
            'TLCA_Concentration': 'TLCA_Conc',
            'GUDCA_Concentration': 'GUDCA_Conc',
            'GCDCA_Concentration': 'GCDCA_Conc',
            'GDCA_Concentration': 'GDCA_Conc',
            'GCA_Concentration': 'GCA_Conc',
            'GLCA_Concentration': 'GLCA_Conc'
        })

        BSH_compound_protein_detailed = BSH_compound_protein_detailed.fillna(0)
        columns_to_round = ['TUDCA', 'TCDCA', 'TDCA', 'TCA', 'TLCA', 'GUDCA', 'GCDCA', 'GDCA', 'GCA', 'GLCA']
        BSH_compound_protein_detailed[columns_to_round] = BSH_compound_protein_detailed[columns_to_round].round(3)

        BSH_compound_protein_detailed = BSH_compound_protein_detailed[BSH_compound_protein_detailed['Replicate'] == 'Average']

        BSH_compound_protein_detailed.to_csv('dataset/generated/BSH_compound_protein_bile_acid_detailed.csv', index=False)

metadata = MetadataGeneration().generate_metadata()
