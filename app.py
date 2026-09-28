#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from flask import Flask, render_template, jsonify

app = Flask(__name__)

TEMPLATE = "TACGGCATTTACGCAACTGATTACAGTCATGCAT"
LENGTH = len(TEMPLATE)
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def generate_simulation_frames():
    """
    Simulates the entire replication pipeline and records a snapshot (frame) 
    at each molecular alteration.
    """
    frames = []
    
    # State tracking structures
    ssbs = set(range(LENGTH))
    RNA_strand = [" "] * LENGTH
    DNA_strand = [" "] * LENGTH
    sealed_nick = False

    def save_frame(enzyme_name, enzyme_marker, enzyme_pos, action_text):
        # Package a deep-copy snapshot of the active molecular environment
        frames.append({
            "template": TEMPLATE,
            "ssbs": ["S" if i in ssbs else " " for i in range(LENGTH)],
            "rna": list(RNA_strand),
            "dna": list(DNA_strand),
            "sealed_nick": sealed_nick,
            "enzyme": {
                "name": enzyme_name,
                "marker": enzyme_marker,
                "pos": enzyme_pos
            },
            "status": action_text
        })

    # Frame 0: Initialization
    save_frame("None", " ", None, "SSBs coated and stabilized single-stranded DNA.")

    fragments_plan = [
        {"name": "Okazaki Fragment 1", "start": 33, "primer_len": 4, "limit": 20},
        {"name": "Okazaki Fragment 2", "start": 19, "primer_len": 4, "limit": 0}
    ]

    # --- PHASE 1: PRIMASE & POLYMERASE III LOOP ---
    for frag in fragments_plan:
        # DNA Primase Entry
        save_frame("DNA Primase", "P", frag['start'], f"DNA Primase binding to initiate {frag['name']}...")
        
        # Build RNA Primer
        for i in range(frag['primer_len']):
            pos = frag['start'] - i
            if pos in ssbs: 
                ssbs.remove(pos)
            base = COMPLEMENT[TEMPLATE[pos]]
            RNA_strand[pos] = 'U' if base == 'T' else base
            save_frame("DNA Primase", "P", pos, "DNA Primase laying down RNA Primer bases (U).")

        # DNA Polymerase III Entry
        pol3_start = frag['start'] - frag['primer_len']
        save_frame("DNA Polymerase III", "3", pol3_start, "DNA Polymerase III docking onto the 3' OH end of the primer.")
        
        # Elongation Phase
        curr_pos = pol3_start
        while curr_pos >= frag['limit']:
            if curr_pos in ssbs: 
                ssbs.remove(curr_pos)
            DNA_strand[curr_pos] = COMPLEMENT[TEMPLATE[curr_pos]]
            save_frame("DNA Polymerase III", "3", curr_pos, f"DNA Pol III extending {frag['name']} with DNA bases.")
            curr_pos -= 1

        save_frame("None", " ", None, f"{frag['name']} synthesis step complete. Notice structural breaks (nicks).")

    # --- PHASE 2: DNA POLYMERASE I MATURATION ---
    for frag in fragments_plan:
        save_frame("DNA Polymerase I", "1", frag['start'], f"DNA Polymerase I targeting RNA primer on {frag['name']}.")
        
        for i in range(frag['primer_len']):
            pos = frag['start'] - i
            rna_base = RNA_strand[pos]
            RNA_strand[pos] = " "  # Excise RNA
            DNA_strand[pos] = 'T' if rna_base == 'U' else rna_base # Replace with DNA
            save_frame("DNA Polymerase I", "1", pos, "DNA Pol I replacing RNA primer base with DNA.")

    # --- PHASE 3: DNA LIGASE BOND REPAIR ---
    nick_position = 20
    save_frame("DNA Ligase", "L", nick_position, "DNA Ligase scanning for structural nicks in the sugar-phosphate backbone...")
    
    sealed_nick = True
    save_frame("DNA Ligase", "L", nick_position, "DNA Ligase catalyzes phosphodiester bond formation! Nick sealed (⦙ removed).")

    # Final wrap frame
    save_frame("None", " ", None, "Ligation complete! The lagging strand is now a single, continuous covalent molecule.")
    
    return frames

@app.route('/')
def index():
    # Serves the main UI page
    return render_template('index.html')

@app.route('/run_simulation')
def run_simulation():
    # Returns the array of states to the browser frontend asynchronously
    simulation_data = generate_simulation_frames()
    return jsonify(simulation_data)

if __name__ == '__main__':
    app.run(debug=True)
