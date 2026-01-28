import streamlit as st
import re
import requests
import pandas as pd
import plotly.graph_objects as go

# --- BIOLOGICAL LOGIC & MOTIF DATABASE ---
DOMAIN_PATTERNS = {
    "KINASE": {
        "pattern": r"G.G..G",  
        "name": "Protein Kinase Motif",
        "function": "Enzymatic Activity",
        "color": "#FF4B4B"
    },
    "SH3": {
        "pattern": r"P..P",    
        "name": "SH3-Binding Motif",
        "function": "Protein Interaction",
        "color": "#1C83E1"
    },
    "WD40": {
        "pattern": r"GH.{20,30}WD", 
        "name": "WD40 Repeat",
        "function": "Scaffolding",
        "color": "#00C04A"
    },
    "N-GLYCO": {
        "pattern": r"N[^P][ST][^P]", 
        "name": "Glycosylation Site",
        "function": "Cell Signaling",
        "color": "#FFAA00"
    }
}

def analyze_sequence(sequence):
    # Clean sequence: remove numbers, spaces, and FASTA headers
    sequence = sequence.upper().strip()
    if sequence.startswith(">"):
        lines = sequence.splitlines()
        sequence = "".join(lines[1:])
    
    sequence = re.sub(r"[\s\d]", "", sequence)

    if not sequence or not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
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

    return {"ptms": ptms, "domains": found_domains, "length": len(sequence)}

def draw_domain_map(domains, seq_length):
    # This initializes the Plotly Figure object
    fig = go.Figure()

    # Draw the main protein backbone
    fig.add_trace(go.Scatter(
        x=[1, seq_length], y=[0, 0],
        mode="lines",
        line=dict(color="#D3D3D3", width=20),
        name="Protein Chain",
        hoverinfo="skip"
    ))

    # Add domains as colored segments
    for d in domains:
        fig.add_trace(go.Bar(
            x=[d["End"] - d["Start"] + 1],
            y=[0],
            base=d["Start"],
            orientation='h',
            marker=dict(color=d["Color"], line=dict(color="black", width=1)),
            name=d["Domain Name"],
            hovertemplate=f"<b>{d['Domain Name']}</b><br>Pos: {d['Start']}-{d['End']}<br>{d['Function']}<extra></extra>"
        ))

    fig.update_layout(
        title="Visual Domain Mapping",
        xaxis=dict(title="Amino Acid Position", range=[0, seq_length + 5], showgrid=True),
        yaxis=dict(showticklabels=False, range=[-1, 1], fixedrange=True),
        height=250,
        margin=dict(l=20, r=20, t=60, b=40),
        showlegend=True,
        plot_bgcolor="white",
        barmode='overlay'
    )
    return fig

# --- STREAMLIT UI ---
st.set_page_config(page_title="Protein Mapper Pro", page_icon="🧬", layout="wide")

st.title("🧬 Protein Domain & PTM Mapper")
st.write("Analyze and visualize protein sequence architecture and post-translational modifications.")

# Initialize session state for input
if 'input_val' not in st.session_state:
    st.session_state.input_val = ""

# Sidebar for example selection
st.sidebar.header("Load Examples")
if st.sidebar.button("Human p53 (Partial)"):
    st.session_state.input_val = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI"
if st.sidebar.button("Src Kinase Segment"):
    st.session_state.input_val = "YVAPSDPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSD"

# Layout for Inputs
col_input, col_uni = st.columns([2, 1])

with col_input:
    seq_input = st.text_area("Amino Acid Sequence:", value=st.session_state.input_val, height=150)

with col_uni:
    uni_id = st.text_input("OR UniProt ID:", placeholder="e.g. P04637")

if st.button("Run Full Analysis", type="primary"):
    final_sequence = seq_input.strip()
    
    # UniProt Fetch Logic
    if uni_id and not final_sequence:
        with st.spinner("Fetching data from UniProt..."):
            try:
                r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uni_id}.fasta")
                if r.status_code == 200:
                    final_sequence = "".join(r.text.split('\n')[1:])
                else:
                    st.error("UniProt ID not found. Please check the ID and try again.")
            except Exception as e:
                st.error(f"Connection Error: {e}")

    if final_sequence:
        res = analyze_sequence(final_sequence)
        
        if "error" in res:
            st.error(res["error"])
        else:
            # 1. Visualization
            st.subheader("📊 Interactive Sequence Map")
            
            if res["domains"]:
                st.plotly_chart(draw_domain_map(res["domains"], res["length"]), use_container_width=True)
            else:
                st.info("No specific functional domains from the database were identified in this sequence.")

            # 2. Result Tables
            tab1, tab2 = st.tabs(["📍 PTM Sites", "🏗️ Structural Motifs"])
            
            with tab1:
                if res["ptms"]:
                    st.dataframe(pd.DataFrame(res["ptms"]), use_container_width=True)
                else:
                    st.info("No common PTM sites (S/T/Y phosphorylation or K acetylation) detected.")
            
            with tab2:
                if res["domains"]:
                    # Remove color column for the display table
                    display_df = pd.DataFrame(res["domains"]).drop(columns=['Color'])
                    st.dataframe(display_df, use_container_width=True)
                else:
                    st.info("No functional motifs identified.")
    else:
        st.warning("Please provide a protein sequence or a UniProt ID.")
