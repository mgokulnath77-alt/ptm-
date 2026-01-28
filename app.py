import streamlit as st
import re
import requests
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# --- BIOLOGICAL LOGIC & MOTIF DATABASE ---
# Added colors for the domain mapping
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

    if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence):
        return {"error": "Invalid characters found. Please use standard amino acid codes."}

    # 1. PTM Identification
    ptms = []
    for i, amino_acid in enumerate(sequence):
        pos = i + 1
        if amino_acid in ['S', 'T', 'Y']:
            ptms.append({"Type": "Phosphorylation", "Residue": amino_acid, "Position": pos})
        if amino_acid == 'K':
            ptms.append({"Type": "Acetylation", "Residue": amino_acid, "Position": pos})

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

    # Draw the main protein backbone
    fig.add_trace(go.Scatter(
        x=[1, seq_length], y=[0, 0],
        mode="lines",
        line=dict(color="#E5ECF6", width=20),
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
            marker=dict(color=d["Color"], line=dict(color="white", width=1)),
            name=d["Domain Name"],
            hovertemplate=f"<b>{d['Domain Name']}</b><br>Pos: {d['Start']}-{d['End']}<br>{d['Function']}<extra></extra>"
        ))

    fig.update_layout(
        title="Interactive Domain Architecture Map",
        xaxis=dict(title="Amino Acid Position", range=[0, seq_length + 10]),
        yaxis=dict(showticklabels=False, range=[-1, 1]),
        height=250,
        margin=dict(l=20, r=20, t=40, b=40),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        plot_bgcolor="white"
    )
    return fig

# --- STREAMLIT UI ---
st.set_page_config(page_title="Protein Profiler Pro", page_icon="🧬", layout="wide")

st.title("🧬 Protein Domain & PTM Mapper")
st.markdown("Visualize the structural architecture and post-translational sites of your protein.")

# Sidebar Examples
st.sidebar.header("Examples")
if st.sidebar.button("Load Human p53"):
    st.session_state.seq = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI"
if st.sidebar.button("Load Src Kinase"):
    st.session_state.seq = "YVAPSDPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSD"

# Input
seq_input = st.text_area("Paste Amino Acid Sequence:", value=st.session_state.get('seq', ''), height=150)
uni_id = st.text_input("OR UniProt ID (e.g. P04637):")

if st.button("Generate Domain Map", type="primary"):
    final_seq = seq_input.strip()
    
    if uni_id and not final_seq:
        with st.spinner("Fetching from UniProt..."):
            r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uni_id}.fasta")
            if r.status_code == 200:
                final_seq = "".join(r.text.split('\n')[1:])
            else:
                st.error("UniProt ID not found.")

    if final_seq:
        res = analyze_sequence(final_seq)
        if "error" in res:
            st.error(res["error"])
        else:
            # 1. Visualization
            st.subheader("📊 Domain Map")
            if res["domains"]:
                st.plotly_chart(draw_domain_map(res["domains"], res["length"]), use_container_width=True)
            else:
                st.warning("No known domains found in this sequence.")

            # 2. Results Tables
            c1, c2 = st.columns(2)
            with c1:
                st.subheader("📍 PTM Sites")
                if res["ptms"]:
                    st.dataframe(pd.DataFrame(res["ptms"]), use_container_width=True)
                else:
                    st.write("No common PTM sites found.")
            
            with c2:
                st.subheader("🏗️ Domain Details")
                if res["domains"]:
                    df_d = pd.DataFrame(res["domains"]).drop(columns=['Color'])
                    st.dataframe(df_d, use_container_width=True)
                else:
                    st.write("No domains identified.")
    else:
        st.info("Please enter a sequence to analyze.")
