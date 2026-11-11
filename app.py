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

# ---------------- EXPANDED COMPOUND DATABASE ----------------

COMPOUNDS = {
    "Water": {
        "smiles": "O",
        "formula": "H2O",
        "uses": "Drinking, cleaning and biological processes."
    },
    "Methane": {
        "smiles": "C",
        "formula": "CH4",
        "uses": "Used as a fuel and natural gas component."
    },
    "Ethanol": {
        "smiles": "CCO",
        "formula": "C2H6O",
        "uses": "Used in sanitizers, solvents and fuels."
    },
    "Glucose": {
        "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "formula": "C6H12O6",
        "uses": "An important energy source for living organisms."
    },
    "Benzene": {
        "smiles": "c1ccccc1",
        "formula": "C6H6",
        "uses": "Used as an industrial chemical starting material."
    },
    "Acetic Acid": {
        "smiles": "CC(=O)O",
        "formula": "C2H4O2",
        "uses": "Main acid present in vinegar."
    },
    "Aspirin": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "formula": "C9H8O4",
        "uses": "A medicine used for pain and fever relief."
    },
    "Acetone": {
        "smiles": "CC(=O)C",
        "formula": "C3H6O",
        "uses": "Common solvent used in nail polish remover and cleaning."
    },
    "Caffeine": {
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "formula": "C8H10N4O2",
        "uses": "Central nervous system stimulant found in coffee and tea."
    },
    "Paracetamol": {
        "smiles": "CC(=O)NC1=CC=C(C=C1)O",
        "formula": "C8H9NO2",
        "uses": "Common pain reliever and fever reducer."
    },
    "Ibuprofen": {
        "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "formula": "C13H18O2",
        "uses": "Nonsteroidal anti-inflammatory drug (NSAID)."
    },
    "Cholesterol": {
        "smiles": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
        "formula": "C27H46O",
        "uses": "Vital lipid component found in cell membranes."
    },
    "Penicillin G": {
        "smiles": "CC1(C(N2C1(C(S2)NC(=O)CC3=CC=CC=C3)C(=O)O)C(=O)O)C",
        "formula": "C16H18N2O4S",
        "uses": "Antibiotic medication used to treat bacterial infections."
    },
    "Urea": {
        "smiles": "NC(=O)N",
        "formula": "CH4N2O",
        "uses": "Used in fertilizers and nitrogen excretion in animals."
    },
    "Formaldehyde": {
        "smiles": "C=O",
        "formula": "CH2O",
        "uses": "Disinfectant, preservative, and chemical precursor."
    },
    "Ammonia": {
        "smiles": "N",
        "formula": "NH3",
        "uses": "Used in agricultural fertilizers and cleaning products."
    },
    "Carbon Dioxide": {
        "smiles": "O=C=O",
        "formula": "CO2",
        "uses": "Greenhouse gas, carbonation in drinks, and fire extinguishers."
    }
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
    - Chemical compound information & SMILES database
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
