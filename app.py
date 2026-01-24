import streamlit as st
import re
import requests
import pandas as pd

# --- BIOLOGICAL LOGIC & MOCK DATABASE ---
DOMAIN_DB = {
    "SH3": {"name": "SH3 Domain", "desc": "Involved in protein-protein interactions.", "function": "Signaling"},
    "KINASE": {"name": "Protein Kinase Domain", "desc": "Catalytic domain for phosphate transfer.", "function": "Enzymatic Activity"},
    "WD40": {"name": "WD40 Repeat", "desc": "Scaffold for multi-protein complex assembly.", "function": "Structural"},
    "HEME": {"name": "Heme-binding Site", "desc": "Binding site for iron-containing porphyrins.", "function": "Metabolism"}
}

def analyze_sequence(sequence):
    sequence = sequence.upper().strip()
    if sequence.startswith(">"):
        sequence = "".join(sequence.splitlines()[1:])
    
    if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
        return {"error": "Invalid characters found in sequence. Please use standard amino acid codes."}

    ptms = []
    for i, amino_acid in enumerate(sequence):
        pos = i + 1
        if amino_acid in ['S', 'T', 'Y']:
            ptms.append({"Type": "Phosphorylation", "Residue": amino_acid, "Position": pos})
        if amino_acid == 'N' and i + 2 < len(sequence):
            if sequence[i+2] in ['S', 'T'] and sequence[i+1] != 'P':
                ptms.append({"Type": "N-glycosylation", "Residue": amino_acid, "Position": pos})
        if amino_acid == 'K':
            ptms.append({"Type": "Acetylation/Ubiquitination", "Residue": amino_acid, "Position": pos})

    found_domains = []
    for motif, info in DOMAIN_DB.items():
        if motif in sequence:
            start = sequence.find(motif) + 1
            found_domains.append({
                "Domain Name": info["name"],
                "Start": start,
                "End": start + len(motif) - 1,
                "Function": info["function"]
            })

    if not found_domains:
        summary = "No known domains identified. Protein might be intrinsically disordered."
    else:
        summary = f"This protein contains {len(found_domains)} identified domains. "
        summary += "Primary functions include: " + ", ".join(set(d['Function'] for d in found_domains)) + "."

    return {"ptms": ptms, "domains": found_domains, "summary": summary}

# --- STREAMLIT INTERFACE ---
st.set_page_config(page_title="Protein Analyzer", page_icon="🧬")

st.title("🧬 Post-Translational Protein Analyzer")
st.markdown("Predict PTMs and identify functional domains from sequences or UniProt IDs.")

# Input Section
seq_input = st.text_area("Input Protein Sequence:", placeholder="e.g., MKWVTFISLLKINASEK")
uni_input = st.text_input("Optional: UniProt ID:", placeholder="e.g., P04637")

if st.button("Analyze Protein"):
    sequence = seq_input.strip()
    
    # UniProt Fetching Logic
    if uni_input and not sequence:
        with st.spinner("Fetching from UniProt..."):
            try:
                resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{uni_input}.fasta")
                if resp.status_code == 200:
                    sequence = "".join(resp.text.split('\n')[1:])
                else:
                    st.error("UniProt ID not found.")
            except:
                st.error("Connection to UniProt failed.")

    if sequence:
        result = analyze_sequence(sequence)
        
        if "error" in result:
            st.error(result["error"])
        else:
            # Display Results
            st.success("Analysis Complete!")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Predicted PTMs")
                if result["ptms"]:
                    st.table(pd.DataFrame(result["ptms"]))
                else:
                    st.info("No PTMs predicted.")

            with col2:
                st.subheader("Identified Domains")
                if result["domains"]:
                    st.table(pd.DataFrame(result["domains"]))
                else:
                    st.info("No domains identified.")

            st.subheader("Functional Summary")
            st.info(result["summary"])
    else:
        st.warning("Please provide a sequence or UniProt ID.")
