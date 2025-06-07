import pandas as pd
import requests

class MetadataGeneration:
    def __init__(self):
        self.bsh_data = None
        self.inhibitors_data = None

    def load_data(self):
        bsh_df = pd.read_csv('dataset/BSH_INHIBITOR_WORKINGDATASET.csv')
        inhibitors_df = pd.read_csv('dataset/inhibitors.csv')
        self.inhibitors_data = inhibitors_df
        self.bsh_data = bsh_df

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

    def get_cid_details(self):
        results = []
        for cid in self.inhibitors_data['CID']:
            core_data = self.fetch_pubchem_compound(cid)
            if core_data:
                results.append(core_data)

        compounds_df = pd.DataFrame(results)
        return compounds_df

    def generate_metadata(self):
        self.load_data()
        compounds = self.get_cid_details()
        compounds_detailed = pd.merge(compounds, self.inhibitors_data, how='left', left_on='cid', right_on='CID')
        BSH_detailed = pd.merge(self.bsh_data, compounds_detailed,  how='left', left_on='Inhibitor', right_on='CODE')
        BSH_detailed['Inhibitor_Input_Concentration'] = 0.1 * BSH_detailed['molecular_weight']
        return BSH_detailed

metadata = MetadataGeneration().generate_metadata()


