import streamlit as st
import pandas as pd
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw
from io import BytesIO
from reportlab.pdfgen import canvas

# ---------------- PAGE CONFIG ----------------

st.set_page_config(
    page_title="ChemiLab Virtual Laboratory",
    page_icon="🧪",
    layout="wide"
)

# ---------------- CUSTOM DESIGN ----------------

st.markdown("""
<style>
.main {
    background-color: #f4f7fb;
}
h1 {
    color: #123c69;
}
.stButton>button {
    background-color: #1769aa;
    color: white;
    border-radius: 8px;
}
</style>
""", unsafe_allow_html=True)

# ---------------- 100 COMPOUNDS DATABASE ----------------

COMPOUNDS = {
    "Water": {"smiles": "O", "formula": "H2O", "uses": "Drinking, cleaning and biological processes."},
    "Methane": {"smiles": "C", "formula": "CH4", "uses": "Used as a fuel and natural gas component."},
    "Ethanol": {"smiles": "CCO", "formula": "C2H6O", "uses": "Used in sanitizers, solvents and fuels."},
    "Glucose": {"smiles": "C(C1C(C(C(C(O1)O)O)O)O)O", "formula": "C6H12O6", "uses": "An important energy source for living organisms."},
    "Benzene": {"smiles": "c1ccccc1", "formula": "C6H6", "uses": "Used as an industrial chemical starting material."},
    "Acetic Acid": {"smiles": "CC(=O)O", "formula": "C2H4O2", "uses": "Main acid present in vinegar."},
    "Aspirin": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "formula": "C9H8O4", "uses": "A medicine used for pain and fever relief."},
    "Acetone": {"smiles": "CC(=O)C", "formula": "C3H6O", "uses": "Common solvent used in nail polish remover."},
    "Caffeine": {"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "formula": "C8H10N4O2", "uses": "Central nervous system stimulant found in coffee and tea."},
    "Paracetamol": {"smiles": "CC(=O)NC1=CC=C(C=C1)O", "formula": "C8H9NO2", "uses": "Common pain reliever and fever reducer."},
    "Ibuprofen": {"smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "formula": "C13H18O2", "uses": "Nonsteroidal anti-inflammatory drug (NSAID)."},
    "Cholesterol": {"smiles": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C", "formula": "C27H46O", "uses": "Vital lipid component found in cell membranes."},
    "Penicillin G": {"smiles": "CC1(C(N2C1(C(S2)NC(=O)CC3=CC=CC=C3)C(=O)O)C(=O)O)C", "formula": "C16H18N2O4S", "uses": "Antibiotic medication used to treat bacterial infections."},
    "Urea": {"smiles": "NC(=O)N", "formula": "CH4N2O", "uses": "Used in fertilizers and nitrogen excretion in animals."},
    "Formaldehyde": {"smiles": "C=O", "formula": "CH2O", "uses": "Disinfectant, preservative, and chemical precursor."},
    "Ammonia": {"smiles": "N", "formula": "NH3", "uses": "Used in agricultural fertilizers and cleaning products."},
    "Carbon Dioxide": {"smiles": "O=C=O", "formula": "CO2", "uses": "Greenhouse gas, carbonation in drinks, and fire extinguishers."},
    "Hydrogen Peroxide": {"smiles": "OO", "formula": "H2O2", "uses": "Antiseptic, bleaching agent, and oxidizing agent."},
    "Sulfuric Acid": {"smiles": "OS(=O)(=O)O", "formula": "H2SO4", "uses": "Used in car batteries, manufacturing, and chemical synthesis."},
    "Hydrochloric Acid": {"smiles": "Cl", "formula": "HCl", "uses": "Stomach acid, metal cleaning, and chemical processing."},
    "Sodium Chloride": {"smiles": "[Na+].[Cl-]", "formula": "NaCl", "uses": "Table salt used for seasoning and preservation."},
    "Sodium Hydroxide": {"smiles": "[Na+].[OH-]", "formula": "NaOH", "uses": "Used in soap making, drain cleaners, and paper production."},
    "Methanol": {"smiles": "CO", "formula": "CH4O", "uses": "Solvent, fuel, and antifreeze component."},
    "Isopropanol": {"smiles": "CC(C)O", "formula": "C3H8O", "uses": "Rubbing alcohol and electronic contact cleaner."},
    "Diethyl Ether": {"smiles": "CCOCC", "formula": "C4H10O", "uses": "Historical anesthetic and laboratory solvent."},
    "Chloroform": {"smiles": "C(Cl)(Cl)Cl", "formula": "CHCl3", "uses": "Historical solvent and reagent in organic chemistry."},
    "Toluene": {"smiles": "Cc1ccccc1", "formula": "C7H8", "uses": "Solvent in paints, thinners, and chemical manufacturing."},
    "Phenol": {"smiles": "Oc1ccccc1", "formula": "C6H6O", "uses": "Used in antiseptics and production of resins/plastics."},
    "Aniline": {"smiles": "Nc1ccccc1", "formula": "C6H7N", "uses": "Precursor to dyes, rubber chemicals, and pharmaceuticals."},
    "Benzoic Acid": {"smiles": "O=C(O)c1ccccc1", "formula": "C7H6O2", "uses": "Food preservative and precursor to chemical synthesis."},
    "Citric Acid": {"smiles": "OC(=O)CC(O)(CC(=O)O)C(=O)O", "formula": "C6H8O7", "uses": "Flavor enhancer and preservative in food and beverages."},
    "Lactic Acid": {"smiles": "CC(O)C(=O)O", "formula": "C3H6O3", "uses": "Found in yogurt, muscle metabolism, and cosmetics."},
    "Ascorbic Acid (Vitamin C)": {"smiles": "C(C1C(C(=C(O1)O)O)=O)O", "formula": "C6H8O6", "uses": "Essential vitamin acting as an antioxidant."},
    "Glycine": {"smiles": "NCC(=O)O", "formula": "C2H5NO2", "uses": "Simplest amino acid found in proteins."},
    "Alanine": {"smiles": "CC(N)C(=O)O", "formula": "C3H7NO2", "uses": "Amino acid used in protein biosynthesis."},
    "Glycerol": {"smiles": "C(C(CO)O)O", "formula": "C3H8O3", "uses": "Humectant, sweetener, and pharmaceutical ingredient."},
    "Propane": {"smiles": "CCC", "formula": "C3H8", "uses": "Fuel for heating and cooking."},
    "Butane": {"smiles": "CCCC", "formula": "C4H10", "uses": "Lighter fluid and portable stove fuel."},
    "Pentane": {"smiles": "CCCCC", "formula": "C5H12", "uses": "Laboratory solvent and blowing agent."},
    "Hexane": {"smiles": "CCCCCC", "formula": "C6H14", "uses": "Extraction solvent for vegetable oils."},
    "Heptane": {"smiles": "CCCCCCC", "formula": "C7H16", "uses": "Solvent and calibration standard for octane rating."},
    "Octane": {"smiles": "CCCCCCCC", "formula": "C8H18", "uses": "Component of gasoline."},
    "Ethylene": {"smiles": "C=C", "formula": "C2H4", "uses": "Plant hormone and polymer manufacturing."},
    "Acetylene": {"smiles": "C#C", "formula": "C2H2", "uses": "Welding fuel and chemical feedstock."},
    "Propylene": {"smiles": "CC=C", "formula": "C3H6", "uses": "Used to manufacture polypropylene plastic."},
    "1-Butene": {"smiles": "CCC=C", "formula": "C4H8", "uses": "Comonomer in polyethylene production."},
    "Isoprene": {"smiles": "CC(=C)C=C", "formula": "C5H8", "uses": "Monomer for synthetic rubber production."},
    "Styrene": {"smiles": "C=Cc1ccccc1", "formula": "C8H8", "uses": "Used to produce polystyrene plastics and resins."},
    "Naphthalene": {"smiles": "c1ccc2ccccc2c1", "formula": "C10H8", "uses": "Mothballs and chemical synthesis precursor."},
    "Anthracene": {"smiles": "c1ccc2cc3ccccc3cc2c1", "formula": "C14H10", "uses": "Dyes and scintillation counter materials."},
    "Pyridine": {"smiles": "c1cc[nH]c1", "formula": "C5H5N", "uses": "Basic heterocyclic organic solvent and reagent."},
    "Pyrrole": {"smiles": "C1=CC=CN1", "formula": "C4H5N", "uses": "Precursor to pharmaceuticals and conducting polymers."},
    "Furan": {"smiles": "C1=CC=CO1", "formula": "C4H4O", "uses": "Solvent and chemical intermediate."},
    "Thiophene": {"smiles": "C1=CC=CS1", "formula": "C4H4S", "uses": "Used in specialty chemical manufacturing."},
    "Imidazole": {"smiles": "C1=CN=CN1", "formula": "C3H4N2", "uses": "Precursor to antifungal medications and ionic liquids."},
    "Oxazole": {"smiles": "C1=NC=CO1", "formula": "C3H3NO", "uses": "Organic synthesis and pharmaceuticals."},
    "Thiazole": {"smiles": "C1=NC=CS1", "formula": "C3H3NS", "uses": "Vitamin B1 synthesis and agrochemicals."},
    "Purine": {"smiles": "c1[nH]c2ncnc2[nH]1", "formula": "C5H4N4", "uses": "Base structure of DNA/RNA nucleobases."},
    "Pyrimidine": {"smiles": "c1cncnc1", "formula": "C4H4N2", "uses": "Found in nucleic acids and pharmaceuticals."},
    "Indole": {"smiles": "c1ccc2c(c1)[nH]cc2", "formula": "C8H7N", "uses": "Fragrances and biochemical precursor (tryptophan)."},
    "Quinoline": {"smiles": "c1ccc2c(c1)nccc2", "formula": "C9H7N", "uses": "Manufacture of antimalarial drugs."},
    "Isoquinoline": {"smiles": "c1ccc2c(c1)cncc2", "formula": "C9H7N", "uses": "Pharmaceutical and alkaloid synthesis."},
    "Acridine": {"smiles": "c1ccc2cc3ccccc3nc2c1", "formula": "C13H9N", "uses": "Dyes and antiseptic agents."},
    "Carbazole": {"smiles": "c1ccc2c(c1)[nH]c3ccccc23", "formula": "C12H9N", "uses": "Photographic plates and pigments."},
    "Dioxane": {"smiles": "C1COCCO1", "formula": "C4H8O2", "uses": "Stable polar aprotic solvent."},
    "Tetrahydrofuran": {"smiles": "C1CCOC1", "formula": "C4H8O", "uses": "Important polar ether solvent (THF)."},
    "Acetonitrile": {"smiles": "CC#N", "formula": "C2H3N", "uses": "Common solvent in HPLC chromatography."},
    "Nitromethane": {"smiles": "C[N+](=O)[O-]", "formula": "CH3NO2", "uses": "Fuel additive and organic solvent."},
    "Carbon Tetrachloride": {"smiles": "C(Cl)(Cl)(Cl)Cl", "formula": "CCl4", "uses": "Historical refrigerant and cleaning solvent."},
    "Dichloromethane": {"smiles": "C(Cl)Cl", "formula": "CH2Cl2", "uses": "Laboratory extraction and paint stripping solvent."},
    "1,2-Dichloroethane": {"smiles": "C(Cl)CCl", "formula": "C2H4Cl2", "uses": "Precursor to vinyl chloride production."},
    "Trichloroethylene": {"smiles": "C(=C(Cl)Cl)Cl", "formula": "C2HCl3", "uses": "Industrial metal degreasing solvent."},
    "Tetrachloroethylene": {"smiles": "C(=C(Cl)Cl)(Cl)Cl", "formula": "C2Cl4", "uses": "Dry cleaning solvent."},
    "Fluorobenzene": {"smiles": "c1ccc(cc1)F", "formula": "C6H5F", "uses": "Pharmaceutical intermediate."},
    "Chlorobenzene": {"smiles": "c1ccc(cc1)Cl", "formula": "C6H5Cl", "uses": "Solvent and intermediate for pesticides/dyes."},
    "Bromobenzene": {"smiles": "c1ccc(cc1)Br", "formula": "C6H5Br", "uses": "Organic synthesis reagent."},
    "Iodobenzene": {"smiles": "c1ccc(cc1)I", "formula": "C6H5I", "uses": "Reagent for aryl-aryl coupling reactions."},
    "Nitrobenzene": {"smiles": "c1ccc(cc1)[N+](=O)[O-]", "formula": "C6H5NO2", "uses": "Production of aniline and lubricants."},
    "Benzaldehyde": {"smiles": "O=Cc1ccccc1", "formula": "C7H6O", "uses": "Almond scent flavoring and chemical precursor."},
    "Acetophenone": {"smiles": "CC(=O)c1ccccc1", "formula": "C8H8O", "uses": "Fragrance and solvent for plastics."},
    "Benzophenone": {"smiles": "O=C(c1ccccc1)c2ccccc2", "formula": "C13H10O", "uses": "UV curing agent and sunscreen ingredient."},
    "Salicylic Acid": {"smiles": "O=C(O)c1ccccc1O", "formula": "C7H6O3", "uses": "Skincare treatments and aspirin precursor."},
    "Phthalic Acid": {"smiles": "O=C(O)c1ccccc1C(=O)O", "formula": "C8H6O4", "uses": "Precursor to plasticizers and polyesters."},
    "Terephthalic Acid": {"smiles": "O=C(O)c1ccc(cc1)C(=O)O", "formula": "C8H6O4", "uses": "Production of PET polyester plastics."},
    "Oxalic Acid": {"smiles": "OC(=O)C(=O)O", "formula": "C2H2O4", "uses": "Rust remover and bleaching agent."},
    "Malonic Acid": {"smiles": "OC(=O)CC(=O)O", "formula": "C3H4O4", "uses": "Organic synthesis building block."},
    "Succinic Acid": {"smiles": "OC(=O)CCC(=O)O", "formula": "C4H6O4", "uses": "Metabolic intermediate and food additive."},
    "Adipic Acid": {"smiles": "OC(=O)CCCCC(=O)O", "formula": "C6H10O4", "uses": "Nylon polymer manufacturing."},
    "Fumaric Acid": {"smiles": "OC(=O)C=CC(=O)O", "formula": "C4H4O4", "uses": "Food acidulant and biochemical intermediate."},
    "Maleic Acid": {"smiles": "OC(=O)C=CC(=O)O", "formula": "C4H4O4", "uses": "Resin and coating manufacturing."},
    "Tartaric Acid": {"smiles": "OC(=O)C(O)C(O)C(=O)O", "formula": "C4H6O6", "uses": "Baking powder component and wine industry."},
    "Valeric Acid": {"smiles": "CCCCC(=O)O", "formula": "C5H10O2", "uses": "Ester production for fragrances."},
    "Caproic Acid": {"smiles": "CCCCCC(=O)O", "formula": "C6H12O2", "uses": "Flavoring agents and chemical synthesis."},
    "Palmitic Acid": {"smiles": "CCCCCCCCCCCCCCCC(=O)O", "formula": "C16H32O2", "uses": "Common saturated fatty acid in animals and plants."},
    "Stearic Acid": {"smiles": "CCCCCCCCCCCCCCCCCC(=O)O", "formula": "C18H36O2", "uses": "Candle making, soaps, and lubricants."},
    "Oleic Acid": {"smiles": "CCCCCCCC\\C=C/CCCCCCCC(=O)O", "formula": "C18H34O2", "uses": "Monounsaturated fatty acid found in olive oil."},
    "Linoleic Acid": {"smiles": "CCCCC\\C=C/C\\C=C/CCCCCCCC(=O)O", "formula": "C18H32O2", "uses": "Essential polyunsaturated omega-6 fatty acid."},
    "Uric Acid": {"smiles": "O=C1NC(=O)C2=C(N1)NC(=N2)NC3=O", "formula": "C5H4N4O3", "uses": "Metabolic breakdown product of purine nucleotides."},
    "Allantoin": {"smiles": "NC(=O)NC1NC(=O)NC1=O", "formula": "C4H6N4O3", "uses": "Skin soothing agent in cosmetics."},
    "Melamine": {"smiles": "NC1=NC(=NC(=N1)N)N", "formula": "C3H6N6", "uses": "Manufacture of laminates, plastics, and glues."}
}

# ---------------- FUNCTIONS ----------------

def get_molecule(name):
    smiles = COMPOUNDS[name]["smiles"]
    return Chem.MolFromSmiles(smiles)


def get_properties(name):
    mol = get_molecule(name)

    return {
        "Compound": name,
        "Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular Weight": round(Descriptors.MolWt(mol), 3),
        "Number of Atoms": mol.GetNumAtoms(),
        "LogP": round(Descriptors.MolLogP(mol), 3),
        "TPSA": round(Descriptors.TPSA(mol), 3),
        "H Bond Donors": Descriptors.NumHDonors(mol),
        "H Bond Acceptors": Descriptors.NumHAcceptors(mol)
    }


def create_pdf(name):
    properties = get_properties(name)
    buffer = BytesIO()

    pdf = canvas.Canvas(buffer)
    pdf.setTitle("ChemiLab Compound Report")

    pdf.setFont("Helvetica-Bold", 18)
    pdf.drawString(50, 790, "ChemiLab - Compound Report")

    pdf.setFont("Helvetica", 12)
    pdf.drawString(50, 750, "Compound: " + name)
    pdf.drawString(50, 725, "Chemical Formula: " + properties["Formula"])
    pdf.drawString(50, 700, "Molecular Weight: " +
                   str(properties["Molecular Weight"]) + " g/mol")
    pdf.drawString(50, 675, "Number of Atoms: " +
                   str(properties["Number of Atoms"]))
    pdf.drawString(50, 650, "LogP: " + str(properties["LogP"]))
    pdf.drawString(50, 625, "TPSA: " + str(properties["TPSA"]))
    pdf.drawString(50, 600, "H Bond Donors: " +
                   str(properties["H Bond Donors"]))
    pdf.drawString(50, 575, "H Bond Acceptors: " +
                   str(properties["H Bond Acceptors"]))

    pdf.drawString(50, 525, "Generated by ChemiLab Virtual Laboratory")

    pdf.save()
    buffer.seek(0)

    return buffer


def detect_functional_groups(mol):
    patterns = {
        "Alcohol": "[OX2H][CX4]",
        "Carboxylic Acid": "[CX3](=O)[OX2H]",
        "Ester": "[CX3](=O)[OX2][#6]",
        "Amine": "[NX3;H2,H1;!$(NC=O)]",
        "Aromatic Ring": "a1aaaaa1",
        "Carbonyl": "[CX3]=[OX1]"
    }

    found = []

    for group, pattern in patterns.items():
        query = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(query):
            found.append(group)

    return found


# ---------------- SIDEBAR ----------------

st.sidebar.title("🧪 ChemiLab")
st.sidebar.caption("Virtual Chemistry Laboratory")

page = st.sidebar.radio(
    "Navigation",
    [
        "Dashboard",
        "Compound Explorer",
        "Molecular Calculator",
        "Structure Viewer",
        "Virtual Experiments",
        "Periodic Table",
        "Chemistry Quiz"
    ]
)

# ---------------- DASHBOARD ----------------

if page == "Dashboard":

    st.title("🧪 ChemiLab")
    st.subheader("Cheminformatics Virtual Laboratory")

    st.write(
        "Explore chemical compounds, calculate molecular properties "
        "and understand molecular structures through an interactive "
        "virtual laboratory."
    )

    col1, col2, col3 = st.columns(3)

    col1.metric("Available Compounds", len(COMPOUNDS))
    col2.metric("Virtual Experiments", 4)
    col3.metric("Periodic Table Elements", "118")

    st.divider()

    st.subheader("What can you explore?")

    st.markdown("""
    - Comprehensive Chemical compound database (100 compounds)
    - Molecular weight & property calculations
    - 2D Molecular structure visualization
    - Functional group identification
    - Compound comparison tools
    - Full Interactive Periodic Table (Elements 1–118)
    - Chemistry knowledge quiz & PDF report generator
    """)

    st.info("Select a module from the sidebar to begin your experiment.")

# ---------------- COMPOUND EXPLORER ----------------

elif page == "Compound Explorer":

    st.title("🔬 Compound Explorer")

    name = st.selectbox("Select a chemical compound", list(COMPOUNDS.keys()))

    mol = get_molecule(name)
    properties = get_properties(name)

    st.subheader(name)

    col1, col2 = st.columns(2)

    with col1:
        st.image(Draw.MolToImage(mol, size=(400, 300)))
    
    with col2:
        st.write("**Chemical Formula:**", properties["Formula"])
        st.write("**Molecular Weight:**", properties["Molecular Weight"], "g/mol")
        st.write("**Number of Atoms:**", properties["Number of Atoms"])
        st.write("**SMILES:**", COMPOUNDS[name]["smiles"])
        st.write("**Common Uses:**", COMPOUNDS[name]["uses"])

    st.subheader("Molecular Properties")
    st.json(properties)

# ---------------- MOLECULAR CALCULATOR ----------------

elif page == "Molecular Calculator":

    st.title("⚗️ Molecular Property Calculator")

    name = st.selectbox("Choose a compound", list(COMPOUNDS.keys()))

    if st.button("Calculate Properties"):

        properties = get_properties(name)

        st.success("Calculation completed successfully!")

        col1, col2, col3 = st.columns(3)

        col1.metric("Molecular Weight", str(properties["Molecular Weight"]) + " g/mol")
        col2.metric("Number of Atoms", properties["Number of Atoms"])
        col3.metric("Molecular Formula", properties["Formula"])

        st.dataframe(
            pd.DataFrame(
                list(properties.items()),
                columns=["Property", "Value"]
            ),
            use_container_width=True
        )

# ---------------- STRUCTURE VIEWER ----------------

elif page == "Structure Viewer":

    st.title("🧬 Molecular Structure Viewer")

    name = st.selectbox("Select Compound", list(COMPOUNDS.keys()))

    mol = get_molecule(name)

    st.subheader(name)

    st.image(
        Draw.MolToImage(mol, size=(600, 400)),
        caption="2D Molecular Structure"
    )

    st.write("**SMILES Representation:**")
    st.code(COMPOUNDS[name]["smiles"])

    st.write("SMILES means Simplified Molecular Input Line Entry System.")

# ---------------- VIRTUAL EXPERIMENTS ----------------

elif page == "Virtual Experiments":

    st.title("🧪 Virtual Chemistry Experiments")

    experiment = st.selectbox(
        "Choose an experiment",
        [
            "Molecular Weight Analysis",
            "Element Composition",
            "Functional Group Identification",
            "Compound Comparison"
        ]
    )

    if experiment == "Molecular Weight Analysis":

        name = st.selectbox("Select Compound", list(COMPOUNDS.keys()))

        if st.button("Run Experiment"):

            properties = get_properties(name)

            st.success("Experiment completed!")

            st.metric(
                "Molecular Weight",
                str(properties["Molecular Weight"]) + " g/mol"
            )

    elif experiment == "Element Composition":

        name = st.selectbox("Select Compound", list(COMPOUNDS.keys()))

        if st.button("Analyze Elements"):

            mol = get_molecule(name)

            elements = {}

            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                elements[symbol] = elements.get(symbol, 0) + 1

            df = pd.DataFrame(
                list(elements.items()),
                columns=["Element", "Number of Atoms"]
            )

            st.dataframe(df, use_container_width=True)

            fig = px.pie(
                df,
                names="Element",
                values="Number of Atoms",
                title="Element Composition"
            )

            st.plotly_chart(fig, use_container_width=True)

    elif experiment == "Functional Group Identification":

        name = st.selectbox("Select Compound", list(COMPOUNDS.keys()))

        if st.button("Identify Groups"):

            mol = get_molecule(name)
            groups = detect_functional_groups(mol)

            if groups:
                st.success("Identified Functional Groups")
                for group in groups:
                    st.write("✓", group)
            else:
                st.info("No supported functional groups detected.")

    elif experiment == "Compound Comparison":

        col1, col2 = st.columns(2)

        with col1:
            first = st.selectbox("First Compound", list(COMPOUNDS.keys()))

        with col2:
            second = st.selectbox(
                "Second Compound",
                list(COMPOUNDS.keys()),
                index=1
            )

        if st.button("Compare Compounds"):

            p1 = get_properties(first)
            p2 = get_properties(second)

            comparison = pd.DataFrame([p1, p2])

            st.dataframe(comparison, use_container_width=True)

            fig = px.bar(
                comparison,
                x="Compound",
                y="Molecular Weight",
                color="Compound",
                title="Molecular Weight Comparison"
            )

            st.plotly_chart(fig, use_container_width=True)

# ---------------- PERIODIC TABLE MODULE (ALL 118 ELEMENTS) ----------------

elif page == "Periodic Table":

    st.title("📊 Interactive Periodic Table (1–118 Elements)")
    st.write(
        "Explore all standard chemical elements using RDKit's built-in periodic table data."
    )

    ptable = Chem.GetPeriodicTable()

    all_elements = []
    for atomic_num in range(1, 119):
        try:
            symbol = ptable.GetElementSymbol(atomic_num)
            weight = ptable.GetAtomicWeight(atomic_num)
            all_elements.append({
                "Atomic Number": atomic_num,
                "Symbol": symbol,
                "Atomic Weight": round(weight, 4)
            })
        except Exception:
            pass

    df_elements = pd.DataFrame(all_elements)

    selected_symbol = st.selectbox(
        "Search or Select Element by Symbol / Atomic Number",
        df_elements["Symbol"] + " (Z=" + df_elements["Atomic Number"].astype(str) + ")"
    )

    selected_num = int(selected_symbol.split("(Z=")[1].replace(")", ""))
    elem_row = df_elements[df_elements["Atomic Number"] == selected_num].iloc[0]

    col1, col2, col3 = st.columns(3)
    col1.metric("Atomic Number", int(elem_row["Atomic Number"]))
    col2.metric("Chemical Symbol", elem_row["Symbol"])
    col3.metric("Atomic Weight", f"{elem_row['Atomic Weight']} u")

    st.subheader("Complete Periodic Elements Database")
    st.dataframe(df_elements, use_container_width=True)

# ---------------- CHEMISTRY QUIZ ----------------

elif page == "Chemistry Quiz":

    st.title("📝 Chemistry Knowledge Quiz")

    questions = [
        {
            "question": "What is the chemical formula of water?",
            "options": ["CO2", "H2O", "O2", "CH4"],
            "answer": "H2O"
        },
        {
            "question": "What is the molecular formula of benzene?",
            "options": ["C6H6", "CH4", "C2H6O", "H2O"],
            "answer": "C6H6"
        },
        {
            "question": "Which compound is commonly found in vinegar?",
            "options": ["Methane", "Benzene", "Acetic Acid", "Glucose"],
            "answer": "Acetic Acid"
        },
        {
            "question": "What does SMILES represent?",
            "options": [
                "A programming language",
                "A molecular representation system",
                "A chemical element",
                "A laboratory instrument"
            ],
            "answer": "A molecular representation system"
        },
        {
            "question": "What is the chemical formula of methane?",
            "options": ["CH4", "CO2", "C6H6", "H2O"],
            "answer": "CH4"
        },
        {
            "question": "Which element is represented by the symbol O?",
            "options": ["Gold", "Oxygen", "Silver", "Iron"],
            "answer": "Oxygen"
        },
        {
            "question": "What is the molecular formula of glucose?",
            "options": ["C6H12O6", "H2O", "CH4", "C6H6"],
            "answer": "C6H12O6"
        },
        {
            "question": "Which software library is used for cheminformatics in this project?",
            "options": ["NumPy", "RDKit", "Matplotlib", "Flask"],
            "answer": "RDKit"
        },
        {
            "question": "What does molecular weight measure?",
            "options": [
                "Molecular mass",
                "Temperature",
                "Pressure",
                "Volume"
            ],
            "answer": "Molecular mass"
        },
        {
            "question": "What is the formula of ethanol?",
            "options": ["C2H6O", "C6H6", "CH4", "CO2"],
            "answer": "C2H6O"
        }
    ]

    with st.form("quiz_form"):

        answers = []

        for i, q in enumerate(questions):

            st.write(f"**Q{i+1}. {q['question']}**")

            answer = st.radio(
                "Choose your answer:",
                q["options"],
                key=f"question_{i}"
            )

            answers.append(answer)

        submitted = st.form_submit_button("Submit Quiz")

    if submitted:

        score = sum(
            answers[i] == questions[i]["answer"]
            for i in range(len(questions))
        )

        st.success(f"Your Score: {score}/10")

        if score >= 8:
            st.balloons()
            st.write("Excellent work!")

        elif score >= 5:
            st.write("Good effort! Keep learning.")

        else:
            st.write("Keep practising and try again.")

        with st.expander("View Correct Answers"):

            for i, q in enumerate(questions):
                st.write(
                    f"Q{i+1}: Correct Answer: {q['answer']}"
                )

# ---------------- PDF REPORT ----------------

st.sidebar.divider()

st.sidebar.subheader("📄 Download Report")

report_compound = st.sidebar.selectbox(
    "Choose a compound for PDF report",
    list(COMPOUNDS.keys())
)

pdf_data = create_pdf(report_compound)

st.sidebar.download_button(
    "Download PDF Report",
    data=pdf_data,
    file_name=f"{report_compound}_Report.pdf",
    mime="application/pdf"
)

st.sidebar.caption("ChemiLab | Academic Virtual lab")
