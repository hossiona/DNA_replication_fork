#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import streamlit as st
import random
from rdkit import Chem
from rdkit.Chem import Draw

# --- Core Data Models ---
class FunctionalGroup:
    def __init__(self, name, formula, smiles, pka, electronegativity, size, 
                 resonance_delocalized=False, hybridization="sp3", inductive_groups=0):
        self.name = name
        self.formula = formula
        self.smiles = smiles  
        self.pka = pka
        self.electronegativity = electronegativity
        self.size = size
        self.resonance_delocalized = resonance_delocalized
        self.hybridization = hybridization
        self.inductive_groups = inductive_groups

# Dataset including Aromatic & Cyclic compounds
DATASET = [
    FunctionalGroup("Ethanol", "CH3CH2OH", "CCO", pka=16.0, electronegativity=3.44, size="small"),
    FunctionalGroup("Acetic Acid", "CH3COOH", "CC(=O)O", pka=4.76, electronegativity=3.44, size="small", resonance_delocalized=True),
    FunctionalGroup("Trichloroacetic Acid", "CCl3COOH", "C(C(=O)O)(Cl)(Cl)Cl", pka=0.66, electronegativity=3.44, size="small", resonance_delocalized=True, inductive_groups=3),
    FunctionalGroup("Ethane", "CH3CH3", "CC", pka=50.0, electronegativity=2.55, size="small", hybridization="sp3"),
    FunctionalGroup("Ethylene", "H2C=CH2", "C=C", pka=44.0, electronegativity=2.55, size="small", hybridization="sp2"),
    FunctionalGroup("Acetylene", "HC#CH", "C#C", pka=25.0, electronegativity=2.55, size="small", hybridization="sp"),
    FunctionalGroup("Ethanethiol", "CH3CH2SH", "CCS", pka=10.6, electronegativity=2.58, size="large"),
    FunctionalGroup("Methylamine", "CH3NH2", "CN", pka=40.0, electronegativity=3.04, size="small"),
    FunctionalGroup("Phenol", "C6H5OH", "C1=CC=C(C=C1)O", pka=10.0, electronegativity=3.44, size="small", resonance_delocalized=True),
    FunctionalGroup("Cyclohexanol", "C6H11OH", "C1CCC(CC1)O", pka=16.0, electronegativity=3.44, size="small", resonance_delocalized=False),
    FunctionalGroup("4-Nitrophenol", "O2NC6H4OH", "C1=CC(=CC=C1[N+](=O)[O-])O", pka=7.15, electronegativity=3.44, size="small", resonance_delocalized=True, inductive_groups=1)
]

def render_rdkit_svg(smiles_string):
    """Generates a native crisp SVG drawing directly parsing the SMILES code."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is not None:
        drawer = Draw.MolDraw2DSVG(280, 200)
        options = drawer.drawOptions()
        options.clearBackground = False  
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    return "<p style='color:red;'>Structure unavailable</p>"

def evaluate_primary_factor(g1, g2):
    """Determines the correct dominant ARIO factor."""
    if g1.size != g2.size:
        winner = g1 if g1.size == "large" else g2
        loser = g2 if winner == g1 else g1
        return winner, "A", f"**Atom Size**: {winner.name} has a larger conjugate base atom (Period 3+), which polarizes and stabilizes negative charge much better than the smaller atom in {loser.name}."
    elif g1.electronegativity != g2.electronegativity:
        winner = g1 if g1.electronegativity > g2.electronegativity else g2
        loser = g2 if winner == g1 else g1
        return winner, "A", f"**Electronegativity**: {winner.name} has a more electronegative atom holding the acidic proton, stabilizing the resulting negative charge better."
    if g1.resonance_delocalized != g2.resonance_delocalized:
        winner = g1 if g1.resonance_delocalized else g2
        return winner, "R", f"**Resonance**: The conjugate base of {winner.name} can delocalize its negative charge across multiple atoms through pi pathways (such as an aromatic benzene ring or carbonyl group), making it vastly more stable."
    hybrid_order = {"sp": 3, "sp2": 2, "sp3": 1}
    if g1.hybridization != g2.hybridization:
        winner = g1 if hybrid_order[g1.hybridization] > hybrid_order[g2.hybridization] else g2
        loser = g2 if winner == g1 else g1
        return winner, "O", f"**Orbital Hybridization**: {winner.name} holds electrons in an {winner.hybridization} orbital. Higher s-character brings the negative charge closer to the nucleus, stabilizing it more than {loser.hybridization}."
    if g1.inductive_groups != g2.inductive_groups:
        winner = g1 if g1.inductive_groups > g2.inductive_groups else g2
        loser = g2 if winner == g1 else g1
        return winner, "I", f"**Inductive Effect**: {winner.name} contains highly electronegative neighboring groups or electron-withdrawing substituents (like a nitro group) that pull electron density away through sigma bonds, dispersing the negative charge."
    winner = g1 if g1.pka < g2.pka else g2
    return winner, "Other", "Subtle structural differences or solvent-stabilization variations."

# --- Streamlit Web App Configuration ---
st.set_page_config(page_title="ARIO Acidity Quiz", page_icon="🧪", layout="centered")

st.title("🧪 The ARIO Acidity Practice Suite")
st.markdown("""
Welcome! This app challenges you to evaluate relative acid strength using the **ARIO** hierarchy:
* **A**tom (Size & Electronegativity) → **R**esonance → **O**rbital → **I**nduction
""")

# Initialize or re-roll the active problem pair
if "g1" not in st.session_state or "g2" not in st.session_state:
    pair = random.sample(DATASET, 2)
    while pair[0].pka == pair[1].pka:
        pair = random.sample(DATASET, 2)
    st.session_state.g1, st.session_state.g2 = pair[0], pair[1]
    st.session_state.submitted = False

g1 = st.session_state.g1
g2 = st.session_state.g2

st.subheader("Compare the Following Pair:")

col1, col2 = st.columns(2)
with col1:
    st.info(f"### Option 1\n**Name:** {g1.name}\n\n**Formula:** `{g1.formula}`")
    st.write(render_rdkit_svg(g1.smiles), unsafe_allow_html=True)
with col2:
    st.success(f"### Option 2\n**Name:** {g2.name}\n\n**Formula:** `{g2.formula}`")
    st.write(render_rdkit_svg(g2.smiles), unsafe_allow_html=True)

st.markdown("---")
student_choice = st.radio("1. Which molecule is the STRONGER acid?", ["Option 1", "Option 2"])
student_factor = st.selectbox(
    "2. Which ARIO factor is the dominant reason for this trend?",
    ["A - Atom (Size or Electronegativity)", "R - Resonance", "O - Orbital Hybridization", "I - Inductive Effect"]
)

true_winner, true_factor_code, explanation = evaluate_primary_factor(g1, g2)
student_picked_group = g1 if student_choice == "Option 1" else g2
student_factor_code = student_factor[0]

if st.button("Submit Answer"):
    st.session_state.submitted = True

if st.session_state.submitted:
    st.markdown("---")
    correct_choice = (student_picked_group == true_winner)
    correct_factor = (student_factor_code == true_factor_code)
    
    if correct_choice and correct_factor:
        st.balloons()
        st.success("🎉 **Perfect!** Both your prediction and chemical reasoning are completely correct.")
    elif correct_choice:
        st.warning("⚠️ **Partial Credit!** You identified the stronger acid, but chose the wrong controlling mechanism.")
    else:
        st.error("❌ **Incorrect.** Review the molecular properties below to see why the stability shifted.")

    st.write("### **Correct Answers Revealed**")
    st.write(f"* **Stronger Acid:** {true_winner.name} (`{true_winner.formula}`)")
    st.write(f"* **Experimental Proof:** pKa of {g1.name} is **{g1.pka}** vs {g2.name} which is **{g2.pka}** (Lower pKa = Stronger Acid)")
    st.write(f"* **Dominant Rule:** {explanation}")
    
    if st.button("Next Problem ➡️"):
        st.session_state.clear()
        st.rerun()
