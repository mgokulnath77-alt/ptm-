from flask import Flask, request, jsonify, render_template_string
import re
import requests

app = Flask(__name__)

# --- BIOLOGICAL LOGIC & MOCK DATABASE ---
# In bioinformatics, domains are conserved parts of a protein sequence 
# that can evolve, function, and exist independently.
DOMAIN_DB = {
    "SH3": {"name": "SH3 Domain", "desc": "Involved in protein-protein interactions.", "function": "Signaling"},
    "KINASE": {"name": "Protein Kinase Domain", "desc": "Catalytic domain for phosphate transfer.", "function": "Enzymatic Activity"},
    "WD40": {"name": "WD40 Repeat", "desc": "Scaffold for multi-protein complex assembly.", "function": "Structural"},
    "HEME": {"name": "Heme-binding Site", "desc": "Binding site for iron-containing porphyrins.", "function": "Metabolism"}
}

def analyze_sequence(sequence):
    """
    Simulates the biological workflow: 
    1. Validation -> 2. PTM Prediction -> 3. Domain ID -> 4. Mapping
    """
    sequence = sequence.upper().strip()
    # Remove FASTA headers if present
    if sequence.startswith(">"):
        sequence = "".join(sequence.splitlines()[1:])
    
    # Basic Protein Validation (Standard 20 Amino Acids)
    if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
        return {"error": "Invalid characters found in sequence. Please use standard amino acid codes."}

    ptms = []
    # 1. PTM Prediction (Rule-based Biochemistry)
    for i, amino_acid in enumerate(sequence):
        pos = i + 1
        
        # Phosphorylation: Ser(S), Thr(T), Tyr(Y) are targets for Kinases
        if amino_acid in ['S', 'T', 'Y']:
            ptms.append({"type": "Phosphorylation", "residue": amino_acid, "pos": pos})
        
        # N-glycosylation: Motif N-X-S/T (where X is not Proline)
        if amino_acid == 'N' and i + 2 < len(sequence):
            if sequence[i+2] in ['S', 'T'] and sequence[i+1] != 'P':
                ptms.append({"type": "N-glycosylation", "residue": amino_acid, "pos": pos})
        
        # Acetylation/Ubiquitination: Lysine (K) is the primary target
        if amino_acid == 'K':
            ptms.append({"type": "Acetylation/Ubiquitination", "residue": amino_acid, "pos": pos})

    # 2. Domain Identification
    found_domains = []
    for motif, info in DOMAIN_DB.items():
        if motif in sequence:
            start = sequence.find(motif) + 1
            found_domains.append({
                "name": info["name"],
                "start": start,
                "end": start + len(motif) - 1,
                "function": info["function"]
            })

    # 3. PTM-Domain Mapping
    # Biological insight: A PTM inside a domain often regulates that domain's function.
    mapping = []
    for p in ptms:
        for d in found_domains:
            if d['start'] <= p['pos'] <= d['end']:
                mapping.append({
                    "ptm": f"{p['type']} ({p['residue']}{p['pos']})",
                    "domain": d['name']
                })

    # 4. Functional Summary
    if not found_domains:
        summary = "No known domains identified. Protein might be intrinsically disordered."
    else:
        summary = f"This protein contains {len(found_domains)} identified domains. "
        summary += "Primary functions include: " + ", ".join(set(d['function'] for d in found_domains)) + "."

    return {
        "sequence": sequence,
        "ptms": ptms,
        "domains": found_domains,
        "mapping": mapping,
        "summary": summary
    }

# --- HTML INTERFACE (Embedded) ---
HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Protein Analyzer</title>
    <style>
        body { font-family: 'Segoe UI', sans-serif; background: #f0f4f8; padding: 40px; color: #333; }
        .card { background: white; max-width: 800px; margin: auto; padding: 30px; border-radius: 12px; box-shadow: 0 4px 15px rgba(0,0,0,0.1); }
        header { border-bottom: 3px solid #27ae60; margin-bottom: 20px; }
        textarea { width: 100%; height: 100px; border: 1px solid #ccc; border-radius: 5px; padding: 10px; box-sizing: border-box; }
        button { background: #27ae60; color: white; border: none; padding: 12px 20px; border-radius: 5px; cursor: pointer; width: 100%; font-size: 16px; margin-top: 10px;}
        .result-box { display: none; margin-top: 25px; border-top: 1px solid #eee; padding-top: 20px; }
        table { width: 100%; border-collapse: collapse; margin: 15px 0; }
        th, td { text-align: left; padding: 10px; border-bottom: 1px solid #eee; }
        th { background: #f9f9f9; color: #2980b9; }
        .summary { background: #e8f4fd; padding: 15px; border-left: 5px solid #2980b9; margin-top: 15px; }
    </style>
</head>
<body>
    <div class="card">
        <header><h1>Post-Translational Protein Analyzer</h1></header>
        <p>Input Protein Sequence (or UniProt ID):</p>
        <textarea id="seqInput" placeholder="Paste sequence here (e.g., MKWVTFISLLKINASEK)"></textarea>
        <input type="text" id="uniInput" placeholder="Optional: UniProt ID (e.g., P04637)" style="width:100%; padding:10px; margin-top:10px; box-sizing:border-box;">
        <button onclick="analyze()">Analyze Protein</button>

        <div id="results" class="result-box">
            <h3>Predicted PTMs</h3>
            <table id="ptmTable"><thead><tr><th>Type</th><th>Site</th></tr></thead><tbody></tbody></table>
            
            <h3>Identified Domains</h3>
            <table id="domTable"><thead><tr><th>Domain</th><th>Range</th><th>Role</th></tr></thead><tbody></tbody></table>
            
            <h3>Functional Summary</h3>
            <div id="summaryText" class="summary"></div>
        </div>
    </div>

    <script>
    async function analyze() {
        const sequence = document.getElementById('seqInput').value;
        const uniprot_id = document.getElementById('uniInput').value;
        
        const response = await fetch('/analyze', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ sequence, uniprot_id })
        });
        
        const data = await response.json();
        if(data.error) return alert(data.error);

        document.getElementById('results').style.display = 'block';
        
        document.querySelector('#ptmTable tbody').innerHTML = data.ptms.map(p => 
            `<tr><td>${p.type}</td><td>${p.residue}${p.pos}</td></tr>`).join('');
            
        document.querySelector('#domTable tbody').innerHTML = data.domains.map(d => 
            `<tr><td>${d.name}</td><td>${d.start}-${d.end}</td><td>${d.function}</td></tr>`).join('');
            
        document.getElementById('summaryText').innerText = data.summary;
    }
    </script>
</body>
</html>
"""

# --- ROUTES ---
@app.route('/')
def home():
    return render_template_string(HTML_TEMPLATE)

@app.route('/analyze', methods=['POST'])
def analyze():
    data = request.get_json()
    sequence = data.get('sequence', '')
    uniprot_id = data.get('uniprot_id', '').strip()

    # External API Call to UniProt if ID is provided
    if uniprot_id and not sequence:
        try:
            resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
            if resp.status_code == 200:
                sequence = "".join(resp.text.split('\\n')[1:])
            else:
                return jsonify({"error": "UniProt ID not found"}), 404
        except:
            return jsonify({"error": "Connection to UniProt failed"}), 500

    if not sequence:
        return jsonify({"error": "Please provide a sequence or ID"}), 400
        
    return jsonify(analyze_sequence(sequence))

if __name__ == '__main__':
    # use_reloader=False prevents the 'signal' error in hosted environments
    app.run(host='0.0.0.0', port=5000, debug=True, use_reloader=False)