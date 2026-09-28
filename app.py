#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from flask import Flask, render_template, jsonify, request
import re

app = Flask(__name__)

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def generate_simulation_frames(template_strand):
    """
    Simulates lagging strand replication with tracking for active Okazaki fragments
    and base pair complement validation.
    """
    frames = []
    length = len(template_strand)
    
    # State tracking structures
    ssbs = set(range(length))
    RNA_strand = [" "] * length
    DNA_strand = [" "] * length
    okazaki_row = [" "] * length  # Track Okazaki fragment IDs visually
    sealed_nick = False

    def save_frame(enzyme_name, enzyme_marker, enzyme_pos, frag_id, action_text):
        frames.append({
            "template": template_strand,
            "ssbs": ["S" if i in ssbs else " " for i in range(length)],
            "rna": list(RNA_strand),
            "dna": list(DNA_strand),
            "okazaki": list(okazaki_row),
            "sealed_nick": sealed_nick,
            "enzyme": {
                "name": enzyme_name,
                "marker": enzyme_marker,
                "pos": enzyme_pos
            },
            "status": action_text
        })

    # Initial snapshot
    save_frame("None", " ", None, " ", "SSBs coated and stabilized the single-stranded template.")

    midpoint = length // 2
    fragments_plan = [
        {"id": "1", "name": "Okazaki Fragment 1", "start": length - 1, "primer_len": min(4, length - midpoint), "limit": midpoint},
        {"id": "2", "name": "Okazaki Fragment 2", "start": midpoint - 1, "primer_len": min(4, midpoint), "limit": 0}
    ]

    # --- PHASE 1: PRIMASE & POLYMERASE III LOOP ---
    for frag in fragments_plan:
        if frag['start'] < 0 or frag['primer_len'] <= 0:
            continue
            
        # DNA Primase Entry
        save_frame("DNA Primase", "P", frag['start'], frag['id'], f"DNA Primase binding to initiate {frag['name']}...")
        
        # Build RNA Primer
        for i in range(frag['primer_len']):
            pos = frag['start'] - i
            if pos < 0: break
            if pos in ssbs: ssbs.remove(pos)
            base = COMPLEMENT[template_strand[pos]]
            RNA_strand[pos] = 'U' if base == 'T' else base
            okazaki_row[pos] = frag['id']  # Mark fragment ownership
            save_frame("DNA Primase", "P", pos, frag['id'], "DNA Primase laying down RNA Primer bases (U).")

        # DNA Polymerase III Entry
        pol3_start = frag['start'] - frag['primer_len']
        if pol3_start >= frag['limit']:
            save_frame("DNA Polymerase III", "3", pol3_start, frag['id'], "DNA Polymerase III docking onto the 3' OH end of the primer.")
        
        # Elongation Phase
        curr_pos = pol3_start
        while curr_pos >= frag['limit']:
            if curr_pos in ssbs: ssbs.remove(curr_pos)
            DNA_strand[curr_pos] = COMPLEMENT[template_strand[curr_pos]]
            okazaki_row[curr_pos] = frag['id']  # Mark fragment ownership
            save_frame("DNA Polymerase III", "3", curr_pos, frag['id'], f"DNA Pol III extending {frag['name']} with DNA bases.")
            curr_pos -= 1

        save_frame("None", " ", None, " ", f"{frag['name']} synthesis step complete. Notice structural breaks (nicks).")

    # --- PHASE 2: DNA POLYMERASE I MATURATION ---
    for frag in fragments_plan:
        if frag['start'] < 0: continue
        save_frame("DNA Polymerase I", "1", frag['start'], frag['id'], f"DNA Polymerase I targeting RNA primer on {frag['name']}.")
        
        for i in range(frag['primer_len']):
            pos = frag['start'] - i
            if pos < 0: break
            rna_base = RNA_strand[pos]
            if rna_base == " ": continue
            RNA_strand[pos] = " "  # Excise RNA
            DNA_strand[pos] = 'T' if rna_base == 'U' else rna_base # Replace with DNA
            save_frame("DNA Polymerase I", "1", pos, frag['id'], "DNA Pol I replacing RNA primer base with DNA.")

    # --- PHASE 3: DNA LIGASE BOND REPAIR ---
    if length > 3:
        nick_position = midpoint
        save_frame("DNA Ligase", "L", nick_position, " ", "DNA Ligase scanning for structural nicks in the sugar-phosphate backbone...")
        
        sealed_nick = True
        save_frame("DNA Ligase", "L", nick_position, " ", "DNA Ligase catalyzes phosphodiester bond formation! Nick sealed (⦙ removed).")

    # Final wrap frame
    save_frame("None", " ", None, " ", "Ligation complete! The lagging strand is now a single, continuous covalent molecule.")
    
    return frames, midpoint

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run_simulation', methods=['POST'])
def run_simulation():
    user_input = request.json.get('template', '').upper().strip()
    
    if not user_input or not re.match("^[ATCG]+$", user_input):
        return jsonify({"error": "Invalid sequence. Please enter character strings containing only A, T, C, and G."}), 400
        
    if len(user_input) < 10 or len(user_input) > 50:
        return jsonify({"error": "For optimal visual presentation, please input a sequence between 10 and 50 bases long."}), 400

    frames, nick_idx = generate_simulation_frames(user_input)
    return jsonify({"frames": frames, "nick_index": nick_idx})

if __name__ == '__main__':
    app.run(debug=True)
