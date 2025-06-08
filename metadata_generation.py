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
        bsh_df = pd.read_csv('dataset/BSH_INHIBITOR_WORKINGDATASET.csv')
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
            "molecular_weight": analysis.molecular_weight(),
            "aromaticity": analysis.aromaticity(),
            "instability_index": analysis.instability_index(),
            "isoelectric_point": analysis.isoelectric_point(),
            "gravy": analysis.gravy(),
            "amino_acid_percent": analysis.get_amino_acids_percent()
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
        compounds_detailed = pd.merge(compounds, self.inhibitors_data, how='left', left_on='cid', right_on='CID')
        BSH_compound_detailed = pd.merge(self.bsh_data, compounds_detailed,  how='left', left_on='Inhibitor', right_on='CODE')
        protein_bsh_detailed = pd.merge(protein_bsh_enrich, self.bsh_enzymes_data, how='left', left_on='accession', right_on='PROTEIN_ACCESSION_NUMBER')
        BSH_compound_protein_detailed = pd.merge(BSH_compound_detailed, protein_bsh_detailed, how='left', left_on='BSH enzyme', right_on='BSH_ENZYME')
        BSH_compound_protein_detailed['Inhibitor_Concentration'] = (
                float(self.micromolar_inhibitor) *
                pd.to_numeric(BSH_compound_protein_detailed['molecular_weight_x'], errors='coerce')
        )
        BSH_compound_protein_detailed['Enzyme_Protein_Concentration'] = (
                float(self.micromolar_enzyme_protein) *
                pd.to_numeric(BSH_compound_protein_detailed['molecular_weight_y'], errors='coerce')
        )
        bile_acids_detailed = pd.merge(bile_acids, self.bile_acid_data, how='left', left_on='cid', right_on='CID')
        bile_acids_detailed['Bile_Acid_Concentration'] = (
                float(self.micromolar_bile_acids) *
                pd.to_numeric(bile_acids_detailed['molecular_weight'], errors='coerce')
        )
        return BSH_compound_protein_detailed

metadata = MetadataGeneration().generate_metadata()


