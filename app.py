import streamlit as st
import re
import requests
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# --- BIOLOGICAL LOGIC & MOTIF DATABASE ---
DOMAIN_PATTERNS = {
    "KINASE": {
        "pattern": r"G.G..G",  
        "name": "Protein Kinase Motif",
        "function": "Enzymatic Activity (Phosphate Transfer)",
        "color": "#FF4B4B"
    },
    "SH3": {
        "pattern": r"P..P",    
        "name": "SH3-Binding Motif",
        "function": "Protein-Protein Interaction",
        "color": "#1C83E1"
    },
    "WD40": {
        "pattern": r"GH.{20,30}WD", 
        "name": "WD40 Repeat",
        "function": "Scaffold / Multi-protein Assembly",
        "color": "#00C04A"
    },
    "N-GLYCO": {
        "pattern": r"N[^P][ST][^P]", 
        "name": "Glycosylation Site",
        "function": "Cell Signaling / Stability",
        "color": "#FFAA00"
    }
}

def analyze_sequence(sequence):
    sequence = sequence.upper().strip()
    if sequence.startswith(">"):
        sequence = "".join(sequence.splitlines()[1:])
    
    # Remove any whitespace or numbers often found in FASTA/GenBank
    sequence = re.sub(r"[\s\d]", "", sequence)

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

    # 2. Domain/Motif Identification
    found_domains = []
    for key, info in DOMAIN_PATTERNS.items():
        matches = re.finditer(info["pattern"], sequence)
        for match in matches:
            found_domains.append({
                "Domain Name": info["name"],
                "Start": match.start() + 1,
                "End": match.end(),
                "Function": info["function"],
                "Color": info["color"]
            })

    # 3. Summary
    if not found_domains:
        summary = "No recognized domains found. The protein may be intrinsically disordered or highly novel."
    else:
        unique_funcs = list(set(d['Function'] for d in found_domains))
        summary = f"Identified {len(found_domains)} structural/functional markers. Roles in: {', '.join(unique_funcs)}."

    return {"ptms": ptms, "domains": found_domains, "summary": summary, "length": len(sequence)}

def plot_domain_map(domains, seq_length):
    """Generates a visual map of the protein domains."""
    fig = go.Figure()

    # Draw the backbone of the protein
    fig.add_trace(go.Scatter(
        x=[1, seq_length], y=[0, 0],
        mode="lines",
        line=dict(color="lightgrey", width=10),
        name="Protein Backbone",
        hoverinfo="skip"
    ))

    # Add each domain as a colored bar
    for d in domains:
        fig.add_trace(go.Bar(
            x=[d["End"] - d["Start"] + 1],
            y=[0.1], # Slightly offset from backbone
            base=d["Start"],
            orientation='h',
            marker=dict(color=d["Color"]),
            name=d["Domain Name"],
            hovertemplate=f"<b>{d['Domain Name']}</b><br>Pos: {d['Start']}-{d['End']}<br>{d['Function']}<extra></extra>"
        ))

    fig.update_layout(
        title="Protein Domain Map (Linear Architecture)",
        xaxis_title="Amino Acid Position",
        yaxis=dict(showticklabels=False, range=[-1, 1]),
        height=300,
        showlegend=True,
        template="plotly_white",
        barmode='overlay'
    )
    return fig

# --- STREAMLIT INTERFACE ---
st.set_page_config(page_title="Protein Profiler Pro", page_icon="🧬", layout="wide")

st.title("🧬 Post-Translational Protein Analyzer")
st.markdown("Analyze protein sequences for PTMs, conserved motifs, and **visual domain mapping**.")

# Sidebar
st.sidebar.header("Test Sequences")
if st.sidebar.button("Load Human p53 (Partial)"):
    st.session_state.seq = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI"
if st.sidebar.button("Load Src Kinase Segment"):
    st.session_state.seq = "YVAPSDPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSD"

# Input
seq_input = st.text_area("Input Protein Sequence:", 
                         value=st.session_state.get('seq', ''),
                         placeholder="Paste amino acids here...")
uni_input = st.text_input("OR Enter UniProt ID:", placeholder="e.g., P04637")

if st.button("Run Full Analysis", type="primary"):
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
            st.success(f"Analysis Complete: Protein Length = {result['length']} aa")
            
            # --- Visual Domain Mapping ---
            st.subheader("🗺️ Sequence Domain Map")
            if result["domains"]:
                fig = plot_domain_map(result["domains"], result["length"])
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No domains to map visually.")

            # --- Data Tables ---
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("📍 Predicted PTM Sites")
                if result["ptms"]:
                    st.dataframe(pd.DataFrame(result["ptms"]), use_container_width=True, height=300)
                else:
                    st.info("No common PTM sites predicted.")

            with col2:
                st.subheader("🏗️ Identified Domains & Motifs")
                if result["domains"]:
                    df_domains = pd.DataFrame(result["domains"]).drop(columns=['Color'])
                    st.dataframe(df_domains, use_container_width=True, height=300)
                else:
                    st.info("No functional domains identified.")

            st.subheader("📝 Functional Summary")
            st.info(result["summary"])
    else:
        st.warning("Please enter a sequence or ID to begin.")
