import streamlit as st
import re
import requests
import pandas as pd

# --- BIOLOGICAL LOGIC & MOTIF DATABASE ---
# Definitions for motifs with associated colors for the map
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
        sequence = "".join(sequence.splitlines()[1:])
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

    return {"ptms": ptms, "domains": found_domains, "length": len(sequence), "seq": sequence}

# --- DOMAIN MAPPING VISUALIZATION ---
def draw_domain_map(domains, seq_length):
    fig = go.Figure()

    # Draw the main protein backbone (the grey line)
    fig.add_trace(go.Scatter(
        x=[1, seq_length], y=[0, 0],
        mode="lines",
        line=dict(color="#D3D3D3", width=15),
        name="Protein Chain",
        hoverinfo="skip"
    ))

    # Add domains as colored segments on top of the backbone
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
        plot_bgcolor="white"
    )
    return fig

# --- STREAMLIT UI ---
st.set_page_config(page_title="Protein Mapper", page_icon="🧬", layout="wide")

st.title("🧬 Protein Domain & PTM Mapper")
st.write("Analyze and visualize protein sequence architecture instantly.")

# Use session state to handle example loading
if 'seq_val' not in st.session_state:
    st.session_state.seq_val = ""

col_in, col_btn = st.columns([4, 1])

with col_btn:
    st.write("### Examples")
    if st.button("Human p53"):
        st.session_state.seq_val = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI"
    if st.button("Src Kinase"):
        st.session_state.seq_val = "YVAPSDPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSD"

with col_in:
    seq_input = st.text_area("Paste Amino Acid Sequence:", value=st.session_state.seq_val, height=120)
    uni_id = st.text_input("OR Enter UniProt ID:")

if st.button("Run Analysis", type="primary"):
    input_data = seq_input.strip()
    
    # Fetch from UniProt if ID is provided
    if uni_id and not input_data:
        with st.spinner("Fetching from UniProt..."):
            try:
                r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uni_id}.fasta")
                if r.status_code == 200:
                    input_data = "".join(r.text.split('\n')[1:])
                else:
                    st.error("UniProt ID not found.")
            except:
                st.error("Connection Error.")

    if input_data:
        res = analyze_sequence(input_data)
        if "error" in res:
            st.error(res["error"])
        else:
            # 1. Visualization
            st.subheader("📊 Interactive Sequence Map")
            
            if res["domains"]:
                st.plotly_chart(draw_domain_map(res["domains"], res["length"]), use_container_width=True)
            else:
                st.info("No conserved motifs from the current database were found.")

            # 2. Results Dataframes
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("📍 Predicted PTM Sites")
                if res["ptms"]:
                    st.dataframe(pd.DataFrame(res["ptms"]), use_container_width=True)
                else:
                    st.write("No common PTM sites detected.")
            
            with c2:
                st.subheader("🏗️ Structural Motifs")
                if res["domains"]:
                    df_d = pd.DataFrame(res["domains"]).drop(columns=['Color'])
                    st.dataframe(df_d, use_container_width=True)
                else:
                    st.write("No domains identified.")
    else:
        st.warning("Please provide a sequence or UniProt ID.")

