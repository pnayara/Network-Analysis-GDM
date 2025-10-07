# 🔬 Functional Network Analysis of Placental Genes in Gestational Diabetes Mellitus



![Python](https://img.shields.io/badge/Python-3.12-blue.svg?logo=python&logoColor=white)
![NetworkX](https://img.shields.io/badge/Bioinformatics-Network%20Analysis-green.svg)
![GSEA](https://img.shields.io/badge/Statistics-Pathway%20Enrichment-orange.svg)


---
## 📖 Summary
This analysis focus into the transcriptomic data from the main study by constructing a **protein-protein interaction (PPI) network** from the **differentially expressed genes (DEGs)** in the **placentas** of patients with **Gestational Diabetes Mellitus** (GDM).


---

## ❓ Problem
What are the characteristics of the protein-protein interaction network formed by differentially expressed genes in the placentas of women with gestational diabetes? 

---

## 🧪 Methodology

**Dataset**: List of DEGs (p < 0.05) from the placental transcriptomics study (GDM vs. Control).

- 🌐 **Network Construction**: The **STRING** database was accessed via its API to retrieve protein-protein interactions with a minimum confidence score of 0.7.  
- 🕸️ **Network Analysis**: The **NetworkX** library was used to construct and visualize the interaction network.  
- 🎯 **Community Detection**: The **Louvain algorithm** (via the *community-louvain* package) was applied to identify densely connected gene modules.  
- 📊 **Hub Identification**: Centrality metrics such as **Degree**, **Betweenness**, and **Closeness** were computed to identify the most influential nodes.  
- 📈 **Enrichment Analysis**: Over-Representation Analysis (ORA) was performed for each community using **Gseapy** against the **Reactome 2022** database to assign biological meaning to the modules.  


---

## 💻 Skills

- 🐍 **Python (v3.12)** → Pandas, NumPy, NetworkX, community-louvain, Gseapy, Matplotlib, Seaborn, Requests  
- 📊 **Statistics** → Network Centrality Metrics (Degree, Betweenness, Closeness), Over-Representation Analysis (ORA)  
- 🔬 **Bioinformatics** → PPI Network Analysis, Community Detection, Pathway Enrichment 

---

## 📊 Results

>This analysis revealed a complex interaction network composed of 818 differentially expressed genes  and 1,064 interactions in the placentas of women with GDM. The network was organized into three main functional communities, revealing the key biological mechanisms affected and the genes that orchestrate these changes.

### Functional Modules

Three primary modules were identified, each associated with a biological pathway critical to GDM pathophysiology:

**Immune Response**: Enriched in pathways for antigen presentation and immune receptor signaling (FcγR, Toll-like), indicating an alteration in the placenta's immune tolerance and inflammatory response.


<p align="center">
  <img src="images\bubbleplot_imune.png"
   width="70%" alt="Enriched Pathways Community 20">
</p>


**Metabolism**: Enriched in pathways for metabolic homeostasis highlighting the metabolic core of the pathology.

<p align="center">
  <img src="images\bubbleplot_Metbolims.png"
   width="70%" alt="Enriched Pathways Community 20">
</p>

**Vascular and Tissue Development**: Enriched in *TGFB/BMP* signaling pathways, which are important for placental growth and vascularization.

<p align="center">
  <img src="images\bubbleplot_vascular.png" width="70%" alt="Enriched Pathways Community 20">
</p>

### Key Genes and Implications

>The analysis identified that the metabolism network plays a key intermediary role, linking the physiological processes of the immune and vascular systems in the placenta.


**Central Hubs: High  Degree and Betweenness Centralities**

- **ENPP1** (Metabolism): Acts as a  bridge in the network, connecting different metabolic pathways. Its function is directly linked to insulin resistance.

- **UBE2N** (Immune): Functions as a signaling node that translates the metabolic stress of GDM into inflammation through the activation of the NF-κB pathway.

- **CDH5** (Vascular): Known as VE-Cadherin, it is essential for vascular integrity. Its position as a hub points to altered angiogenesis and vascular function in the GDM placenta.


<p align="center">
  <img src="images\network_degree.png"
   width="70%" alt="Enriched Pathways Community 20">
  <img src="images\network_between.png" width="70%" alt="Enriched Pathways Community 20">
</p>

**Hub with High LogFC**
- **BMP7** (Vascular): being both a hub and a highly expressed gene. This suggests it is a key factor that actively drives pathological structural and vascular changes in the GDM placenta.

<p align="center">
  <img src="images\network_FC.png" width="70%" alt="Enriched Pathways Community 20">
</p>

**Highly LogFC Genes** 
- **IDO2** (Metabolism): Its high expression suggests a strong dysregulation of the tryptophan pathway, impacting the placental inflammatory environment.
- **CARD14** (Immune): A known activator of the pro-inflammatory NF-κB pathway. Its increased expression indicates an amplification of the inflammatory response.

<p align="center">
  <img src="images\scatter_FC_DG Imune.png"
   width="30%" alt="Enriched Pathways Community 20">
   <img src="images\scatter_FC_DG metabolims.png" width="30%" alt="Enriched Pathways Community 20">
   <img src="images\scatter_FC_DG_vascular.png" width="30%" alt="Enriched Pathways Community 20">
</p>

---

## 📝 Conclusion
- The GDM placenta exhibits a pathological interplay among metabolic, immune, and vascular pathways, with metabolism acting as the central mediator.

- The results suggest that metabolic stress triggers an inflammatory state, disrupting vascular development and contributing to GDM-associated placental dysfunction. 

