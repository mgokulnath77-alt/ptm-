import streamlit as st
import re
import requests
import pandas as pd

# --- BIOLOGICAL LOGIC & MOTIF DATABASE ---
# Updated to use Regex patterns for real biological domain detection
DOMAIN_PATTERNS = {
    "KINASE": {
        "pattern": r"G.G..G",  # Simplified Glycine-rich loop (P-loop) common in Kinases
        "name": "Protein Kinase Motif",
        "function": "Enzymatic Activity (Phosphate Transfer)"
    },
    "SH3": {
        "pattern": r"P..P",    # Proline-rich motif often recognized by SH3 domains
        "name": "SH3-Binding Motif",
        "function": "Protein-Protein Interaction"
    },
    "WD40": {
        "pattern": r"GH.{20,30}WD", # Classic GH...WD repeat structure
        "name": "WD40 Repeat",
        "function": "Scaffold / Multi-protein Assembly"
    },
    "N-GLYCO": {
        "pattern": r"N[^P][ST][^P]", # N-glycosylation consensus: N-X-S/T
        "name": "Glycosylation Site",
        "function": "Cell Signaling / Stability"
    }
}

def analyze_sequence(sequence):
    sequence = sequence.upper().strip()
    if sequence.startswith(">"):
        sequence = "".join(sequence.splitlines()[1:])
    
    if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
        return {"error": "Invalid characters found. Please use standard amino acid codes."}

    # 1. PTM Identification
    ptms = []
    for i, amino_acid in enumerate(sequence):
        pos = i + 1
        if amino_acid in ['S', 'T', 'Y']:
            ptms.append({"Type": "Phosphorylation", "Residue": amino_acid, "Position": pos})
        if amino_acid == 'K':
            ptms.append({"Type": "Acetylation/Ubiquitination", "Residue": amino_acid, "Position": pos})

    # 2. Domain/Motif Identification (Using Regex Search)
    found_domains = []
    for key, info in DOMAIN_PATTERNS.items():
        matches = re.finditer(info["pattern"], sequence)
        for match in matches:
            found_domains.append({
                "Domain Name": info["name"],
                "Start": match.start() + 1,
                "End": match.end(),
                "Function": info["function"]
            })

    # 3. Functional Summary Logic
    if not found_domains:
        summary = "No recognized domains found. The protein may be intrinsically disordered or highly novel."
    else:
        # Collect unique functions
        unique_funcs = list(set(d['Function'] for d in found_domains))
        summary = f"Identified {len(found_domains)} structural/functional markers. "
        summary += f"This suggests roles in: {', '.join(unique_funcs)}."

    return {"ptms": ptms, "domains": found_domains, "summary": summary}

# --- STREAMLIT INTERFACE ---
st.set_page_config(page_title="Protein Profiler", page_icon="🧬", layout="wide")

st.title("🧬 Post-Translational Protein Analyzer")
st.markdown("Analyze protein sequences for PTMs and conserved functional motifs.")

# Sidebar for example data
st.sidebar.header("Test Sequences")
if st.sidebar.button("Load Human p53 (Partial)"):
    st.session_state.seq = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI"
if st.sidebar.button("Load Src Kinase Segment"):
    st.session_state.seq = "YVAPSDPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSD"

# Input Section
seq_input = st.text_area("Input Protein Sequence:", 
                         value=st.session_state.get('seq', ''),
                         placeholder="Paste amino acids here...")
uni_input = st.text_input("OR Enter UniProt ID:", placeholder="e.g., P04637")

if st.button("Run Full Analysis"):
    sequence = seq_input.strip()
    
    if uni_input and not sequence:
        with st.spinner("Accessing UniProt Database..."):
            try:
                resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{uni_input}.fasta")
                if resp.status_code == 200:
                    sequence = "".join(resp.text.split('\n')[1:])
                else:
                    st.error("UniProt ID not found.")
            except:
                st.error("UniProt API connection failed.")

    if sequence:
        result = analyze_sequence(sequence)
        
        if "error" in result:
            st.error(result["error"])
        else:
            st.success("Analysis Successfully Executed")
            
            # Use columns for layout
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("📍 Predicted PTM Sites")
                if result["ptms"]:
                    st.dataframe(pd.DataFrame(result["ptms"]), use_container_width=True)
                else:
                    st.info("No common PTM sites predicted.")

            with col2:
                st.subheader("🏗️ Identified Domains & Motifs")
                if result["domains"]:
                    st.dataframe(pd.DataFrame(result["domains"]), use_container_width=True)
                else:
                    st.info("No functional domains identified.")

            st.subheader("📝 Functional Summary")
            st.info(result["summary"])
    else:
        st.warning("Please enter a sequence or ID to begin.")
