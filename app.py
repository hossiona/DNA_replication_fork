#!/usr/bin/env python
# coding: utf-8

# In[ ]:

 
import os
import random
from flask import Flask, render_template, request, session, redirect, url_for

app = Flask(__name__)
# Secure fallback for cookie session validation
app.secret_key = os.environ.get("FLASK_SECRET_KEY", "ario_teaching_secret_key_12345")

class FunctionalGroup:
    def __init__(self, id_str, name, formula, pka, electronegativity, size, 
                 resonance_delocalized=False, hybridization="sp3", inductive_groups=0):
        self.id = id_str
        self.name = name
        self.formula = formula
        self.pka = pka  
        self.electronegativity = electronegativity  
        self.size = size                            
        self.resonance_delocalized = resonance_delocalized  
        self.hybridization = hybridization          
        self.inductive_groups = inductive_groups    

    def to_dict(self):
        return self.__dict__

# Static Curated Dataset Map
DATASET = {
    "ethanol": FunctionalGroup("ethanol", "Ethanol", "CH3CH2OH", 16.0, 3.44, "small"),
    "acetic_acid": FunctionalGroup("acetic_acid", "Acetic Acid", "CH3COOH", 4.76, 3.44, "small", resonance_delocalized=True),
    "tca": FunctionalGroup("tca", "Trichloroacetic Acid", "CCl3COOH", 0.66, 3.44, "small", resonance_delocalized=True, inductive_groups=3),
    "ethane": FunctionalGroup("ethane", "Ethane", "CH3CH3", 50.0, 2.55, "small", hybridization="sp3"),
    "ethylene": FunctionalGroup("ethylene", "Ethylene", "H2C=CH2", 44.0, 2.55, "small", hybridization="sp2"),
    "acetylene": FunctionalGroup("acetylene", "Acetylene", "HC#CH", 25.0, 2.55, "small", hybridization="sp"),
    "thiol": FunctionalGroup("thiol", "Ethanethiol", "CH3CH2SH", 10.6, 2.58, "large"),
    "methylamine": FunctionalGroup("methylamine", "Methylamine", "CH3NH2", 40.0, 3.04, "small")
}

def evaluate_primary_factor(g1, g2):
    """Determines the correct dominant ARIO factor."""
    if g1.size != g2.size:
        winner = g1 if g1.size == "large" else g2
        loser = g2 if winner == g1 else g1
        return winner.id, "A", f"**Atom Size**: {winner.name} has a larger conjugate base atom (Period 3+), which polarizes and stabilizes negative charge much better than the smaller atom in {loser.name}."
    elif g1.electronegativity != g2.electronegativity:
        winner = g1 if g1.electronegativity > g2.electronegativity else g2
        loser = g2 if winner == g1 else g1
        return winner.id, "A", f"**Electronegativity**: {winner.name} has a more electronegative atom holding the acidic proton, stabilizing the resulting negative charge better."

    if g1.resonance_delocalized != g2.resonance_delocalized:
        winner = g1 if g1.resonance_delocalized else g2
        return winner.id, "R", f"**Resonance**: The conjugate base of {winner.name} can delocalize its negative charge across multiple atoms through pi pathways, making it vastly more stable."

    hybrid_order = {"sp": 3, "sp2": 2, "sp3": 1}
    if g1.hybridization != g2.hybridization:
        winner = g1 if hybrid_order[g1.hybridization] > hybrid_order[g2.hybridization] else g2
        loser = g2 if winner == g1 else g1
        return winner.id, "O", f"**Orbital Hybridization**: {winner.name} holds electrons in an {winner.hybridization} orbital. Higher s-character brings the negative charge closer to the nucleus, stabilizing it more than {loser.hybridization}."

    if g1.inductive_groups != g2.inductive_groups:
        winner = g1 if g1.inductive_groups > g2.inductive_groups else g2
        loser = g2 if winner == g1 else g1
        return winner.id, "I", f"**Inductive Effect**: {winner.name} contains highly electronegative neighboring groups that pull electron density away through sigma bonds, dispersing the negative charge."

    winner = g1 if g1.pka < g2.pka else g2
    return winner.id, "Other", "Subtle structural differences or solvent-stabilization variations."

@app.route("/", methods=["GET", "POST"])
def index():
    # Setup initial session pair if empty
    if "g1_id" not in session or "g2_id" not in session:
        keys = list(DATASET.keys())
        p1, p2 = random.sample(keys, 2)
        while DATASET[p1].pka == DATASET[p2].pka:
            p1, p2 = random.sample(keys, 2)
        session["g1_id"] = p1
        session["g2_id"] = p2
        session["submitted"] = False

    g1 = DATASET[session["g1_id"]]
    g2 = DATASET[session["g2_id"]]
    
    true_winner_id, true_factor_code, explanation = evaluate_primary_factor(g1, g2)
    
    feedback = None
    
    if request.method == "POST" and "submit_btn" in request.form:
        session["submitted"] = True
        student_choice_id = request.form.get("student_choice")
        student_factor_code = request.form.get("student_factor")
        
        correct_choice = (student_choice_id == true_winner_id)
        correct_factor = (student_factor_code == true_factor_code)
        
        if correct_choice and correct_factor:
            feedback = {"status": "success", "msg": "🎉 Perfect! Both your prediction and chemical reasoning are completely correct."}
        elif correct_choice:
            feedback = {"status": "warning", "msg": "⚠️ Partial Credit! You identified the stronger acid, but chose the wrong controlling mechanism."}
        else:
            feedback = {"status": "danger", "msg": "❌ Incorrect. Review the molecular properties below to see why the stability shifted."}

    return render_template(
        "index.html", 
        g1=g1, g2=g2, 
        submitted=session.get("submitted", False),
        feedback=feedback,
        true_winner=DATASET[true_winner_id],
        explanation=explanation
    )

@app.route("/next")
def next_question():
    session.pop("g1_id", None)
    session.pop("g2_id", None)
    session["submitted"] = False
    return redirect(url_for("index"))

if __name__ == "__main__":
    app.run(debug=True)
