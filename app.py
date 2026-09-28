#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from flask import Flask, render_template, jsonify, request
import re
import random

app = Flask(__name__)

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def generate_full_simulation(template_strand, trigger_mutation=True):
    """
    Simulates full replication with Helicase unzipping, SSB binding, 
    Pol III misincorporation & proofreading, Okazaki tracking, and Ligase action.
    """
    frames = []
    length = len(template_strand)
    
    # State structures
    unzipped_mask = [False] * length  # True means Helicase broke the bonds here
    ssbs = set()
    RNA_strand = [" "] * length
    DNA_strand = [" "] * length
    okazaki_row = [" "] * length
    sealed_nick = False
    
    # Energy Counters
    energy = {"ATP": 0, "GTP": 0, "dNTP": 0}
    
    # Midpoint split for 2 Okazaki fragments
    midpoint = length // 2

    def save_frame(enzyme_name, enzyme_marker, enzyme_pos, action_text):
        frames.append({
            "template": template_strand,
            "unzipped": list(unzipped_mask),
            "ssbs": ["S" if i in ssbs else " " for i in range(length)],
            "rna": list(RNA_strand),
            "dna": list(DNA_strand),
            "okazaki": list(okazaki_row),
            "sealed_nick": sealed_nick,
            "energy": dict(energy),
            "enzyme": {"name": enzyme_name, "marker": enzyme_marker, "pos": enzyme_pos},
            "status": action_text
        })

    # --- PHASE 1: HELICASE & SSB ACTIVATION ---
    # Helicase unzips the downstream portion first so Fragment 1 can start
    save_frame("None", " ", None, "Double-stranded DNA genomic duplex ready for synthesis.")
    
    for pos in range(midpoint, length):
        unzipped_mask[pos] = True
        ssbs.add(pos)
        energy["ATP"] += 1  # 1 ATP per base unzipped
        save_frame("DNA Helicase", "H", pos, f"Helicase consuming ATP to split hydrogen bonds at position {pos}.")

    # --- PHASE 2: SYNTHESIZING FRAGMENT 1 ---
    frag1_start = length - 1
    frag1_primer_len = min(4, length - midpoint)
    frag1_limit = midpoint

    # Primase
    save_frame("DNA Primase", "P", frag1_start, "DNA Primase binding to downstream template initiating Fragment 1...")
    for i in range(frag1_primer_len):
        pos = frag1_start - i
        if pos in ssbs: ssbs.remove(pos)
        base = COMPLEMENT[template_strand[pos]]
        RNA_strand[pos] = 'U' if base == 'T' else base
        okazaki_row[pos] = '1'
        energy["GTP"] += 1  # RNA synthesis consumes GTP/NTPs
        save_frame("DNA Primase", "P", pos, "Primase deploying RNA Primer nucleotides.")

    # Pol III Elongation (Normal)
    pol3_pos = frag1_start - frag1_primer_len
    save_frame("DNA Polymerase III", "3", pol3_pos, "DNA Polymerase III docking onto Fragment 1 primer.")
    while pol3_pos >= frag1_limit:
        if pol3_pos in ssbs: ssbs.remove(pos)
        DNA_strand[pol3_pos] = COMPLEMENT[template_strand[pol3_pos]]
        okazaki_row[pol3_pos] = '1'
        energy["dNTP"] += 1
        save_frame("DNA Polymerase III", "3", pol3_pos, "DNA Pol III expanding Fragment 1.")
        pol3_pos -= 1

    # --- PHASE 3: HELICASE UNZIPS REMAINDER FOR FRAGMENT 2 ---
    save_frame("None", " ", None, "Replication fork moving forward. Unzipping rest of the sequence...")
    for pos in range(0, midpoint):
        unzipped_mask[pos] = True
        ssbs.add(pos)
        energy["ATP"] += 1
        save_frame("DNA Helicase", "H", pos, f"Helicase unzipping upstream sequence at position {pos}.")

    # --- PHASE 4: SYNTHESIZING FRAGMENT 2 + PROOFREADING LOOP ---
    frag2_start = midpoint - 1
    frag2_primer_len = min(4, midpoint)
    frag2_limit = 0

    # Primase Fragment 2
    for i in range(frag2_primer_len):
        pos = frag2_start - i
        if pos in ssbs: ssbs.remove(pos)
        base = COMPLEMENT[template_strand[pos]]
        RNA_strand[pos] = 'U' if base == 'T' else base
        okazaki_row[pos] = '2'
        energy["GTP"] += 1
        save_frame("DNA Primase", "P", pos, "Primase laying down RNA primer for Fragment 2.")

    # Pol III Elongation with Mutation/Proofreading trigger
    pol3_pos = frag2_start - frag2_primer_len
    
    # We choose a specific spot in fragment 2 elongation to trigger a mutation frame
    mutation_spot = frag2_limit + (pol3_start := pol3_pos - frag2_limit) // 2
    
    while pol3_pos >= frag2_limit:
        if pol3_pos in ssbs: ssbs.remove(pol3_pos)
        
        if trigger_mutation and pol3_pos == mutation_spot:
            # Inject a deliberate mismatched base (e.g. force an 'A' regardless of complement rules)
            correct_base = COMPLEMENT[template_strand[pol3_pos]]
            wrong_base = 'A' if correct_base != 'A' else 'C'
            
            DNA_strand[pol3_pos] = wrong_base
            okazaki_row[pol3_pos] = '2'
            energy["dNTP"] += 1
            save_frame("DNA Polymerase III", "3", pol3_pos, "⚠️ WARNING: DNA Pol III structural slip! Mismatched base incorporated.")
            
            # Proofreading step activation frame
            save_frame("DNA Polymerase III (Exonuclease)", "E", pol3_pos, "❌ Mismatch detected by confirmation check! Activating 3'→5' Exonuclease...")
            
            # Remove mismatch
            DNA_strand[pol3_pos] = " "
            okazaki_row[pol3_pos] = " "
            save_frame("DNA Polymerase III (Exonuclease)", "E", pol3_pos, "Exonuclease excised the erroneous nucleotide base.")
            
            # Re-synthesize correctly
            DNA_strand[pol3_pos] = correct_base
            okazaki_row[pol3_pos] = '2'
            energy["dNTP"] += 1
            save_frame("DNA Polymerase III", "3", pol3_pos, "Pol III re-synthesized the position properly via high-fidelity match.")
        else:
            DNA_strand[pol3_pos] = COMPLEMENT[template_strand[pol3_pos]]
            okazaki_row[pol3_pos] = '2'
            energy["dNTP"] += 1
            save_frame("DNA Polymerase III", "3", pol3_pos, "DNA Pol III extending Fragment 2.")
            
        pol3_pos -= 1

    # --- PHASE 5: POL I MATURATION ---
    plans = [{"id": '1', "start": frag1_start, "len": frag1_primer_len}, {"id": '2', "start": frag2_start, "len": frag2_primer_len}]
    for p in plans:
        save_frame("DNA Polymerase I", "1", p['start'], f"DNA Pol I binding to remove RNA primer on Fragment {p['id']}.")
        for i in range(p['len']):
            pos = p['start'] - i
            rna_b = RNA_strand[pos]
            if rna_b == " ": continue
            RNA_strand[pos] = " "
            DNA_strand[pos] = 'T' if rna_b == 'U' else rna_b
            save_frame("DNA Polymerase I", "1", pos, "DNA Pol I replacing RNA primer nucleotide with DNA base.")

    # --- PHASE 6: LIGASE BACKBONE SEAL ---
    if length > 3:
        nick_pos = midpoint
        save_frame("DNA Ligase", "L", nick_pos, "DNA Ligase targeting structural nick between Fragment 2 and Fragment 1.")
        sealed_nick = True
        energy["ATP"] += 2  # Ligase consumes ATP to drive phosphodiester condensation reactions
        save_frame("DNA Ligase", "L", nick_pos, "DNA Ligase consumed ATP. Phosphodiester backbone completed! Nick sealed.")

    save_frame("None", " ", None, "Replication process finalized successfully.")
    return frames, midpoint

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run_simulation', methods=['POST'])
def run_simulation():
    user_input = request.json.get('template', '').upper().strip()
    mutate_opt = request.json.get('mutate', True)
    
    if not user_input or not re.match("^[ATCG]+$", user_input):
        return jsonify({"error": "Invalid sequence. Enter character strings containing only A, T, C, and G."}), 400
    if len(user_input) < 12 or len(user_input) > 50:
        return jsonify({"error": "Please enter a DNA template between 12 and 50 base pairs for layout sizing."}), 400

    frames, nick_idx = generate_full_simulation(user_input, mutate_opt)
    return jsonify({"frames": frames, "nick_index": nick_idx})

if __name__ == '__main__':
    app.run(debug=True)
