# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 11:58:43 2025

@author: Nayara
"""
#%%Install
#!pip install pyvis

#%%Packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import requests # HTTP  API  STRING
import io # API to file
import networkx as nx # netwok 
import community as community_louvain # communities
import gseapy as gp # Enrichment


#%% Data loading and data wrangling
degs = pd.read_csv('DEGs_placenta_limmaVoon.csv')
degs.head()

# Filter p-value
degs_sig =degs[degs['P.Value'] < 0.05]
print(degs_sig['P.Value'].head())
print(degs_sig['P.Value'].tail())
len(degs_sig)
degs_sig['Symbol']

#%% STRIG API 

# Add '\r' to STRING API parse

degs_sig['Symbol']
gene_symbol = degs_sig['Symbol'].to_list()
gene_string = '\r'.join(gene_symbol)
gene_string

#String API
string_url_base = 'https://string-db.org/api'
output_format = 'tsv'
method = 'network'

request_url = f"{string_url_base}/{output_format}/{method}"
print(f"URL da requisição: {request_url}")

# Parameters for requests
params ={
    'identifiers': gene_string,
    'species': 9606,
    'required_score' : 400,
    'caller_identity': 'String_DEGs_Placenta'
}

# Request API
response = requests.post(request_url, data=params, timeout=60)
string_data = response.text
response.status_code

#%% Data from String request
string_df =pd.read_csv(io.StringIO(string_data), sep ='\t')
string_df.sample(10)

# Scores
string_df.sort_values(by='score', ascending=False)

string_df.columns

#filter for score/ confidence minimum 0.7

string_df_filt = string_df[string_df['score']>= 0.7]
string_df_filt = string_df_filt[['preferredName_A', 'preferredName_B', 'score','escore']]
string_df_filt.head()

string_df_filt.info()

#%% NETWORK analysis

#nodes and edges
source = string_df_filt['preferredName_A']
target = string_df_filt['preferredName_B']
weights = string_df_filt['score']

#checking
edge_data = list(zip(source, target, weights))
edge_list = edge_data[:10]
print(edge_list)


# Create graph
G = nx.Graph()
for src, tgt, weight in edge_data:
    G.add_edge(src, tgt, weight=weight)


plt.figure(figsize=(12, 8))
pos =  nx.kamada_kawai_layout( G, 
                              weight='weight', 
                              scale=4, 
                              center=(0, 0))
nx.draw(G, pos, 
        with_labels=False, 
        node_color='lightblue', 
        node_size=100, 
        font_size=8,
        font_weight='bold')
plt.title("")
plt.show()


#%% Communities 

G= nx.Graph()
for src, tgt, w in edge_data:
    G.add_edge(src, tgt, weight = w) 

#Detecting communities
partition = community_louvain.best_partition(G, random_state=42)
nx.set_node_attributes(G, partition , 'community') 

#Id hubs -centrality
degree_centrality = nx.degree_centrality(G)

#betweenness Centrality
betweenness_centrality = nx.betweenness_centrality(G, normalized=True)

#Closenesse Centrality 
closeness_centrality = nx.closeness_centrality(G)

# df communities and centrality
communities_df = pd.DataFrame({'community':partition,
                        'degree_centrality':degree_centrality,
                        'betweenness':betweenness_centrality,
                        'closeness':closeness_centrality
                        })
communities_df.reset_index(inplace = True)
communities_df['community']= communities_df['community'].astype('str')

#communities sizes
communities_size = communities_df['community'].value_counts()
communities_size = communities_size.to_frame()
communities_size.reset_index(inplace = True)

# Filter for communities 
commu_filter = communities_size[
    (communities_size['count'] >=10) &
     (communities_size['count'] <= 60)]

#Genes, communities and Centrality 
community_filter = communities_df[communities_df['community'].isin(commu_filter['community'])]
community_filter.rename(columns = {'index':'genes'}, inplace = True)

#network of community genes
mask = (string_df_filt['preferredName_A'].isin(community_filter['genes']) |
        string_df_filt['preferredName_B'].isin(community_filter['genes']))

network_community = string_df_filt[mask]    

#%% Pathway analysis for community genes 
        
#dict for genes in each community

dict_comm = community_filter.groupby('community')['genes'].agg(list).to_dict()

#datasets for enrichment 
gp.get_library_name(organism='Human')
#gene_sets:'GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022'

#ORA analysis
ora = []
for comm, genes in dict_comm.items():
    enr = gp.enrich(gene_list= genes,
                    gene_sets='Reactome_2022',
                    outdir = None)
    results_ora = enr.results
    results_ora['Community'] = comm
    ora.append(results_ora)
    
ora_communities = pd.concat(ora, ignore_index=True)

#Selecting Enrich Adjusted-value <0.05
ora_communities_filt = ora_communities[ora_communities['Adjusted P-value']<=0.05]

#Group communities pathways
df_ora_comm = ora_communities_filt.groupby('Community')['Term'].agg(list)
df_ora_comm = df_ora_comm.to_frame()
df_ora_comm.reset_index(inplace = True)
#df_ora_comm.to_excel('Communities_enrichPathways.xlsx', index = False)


# PLOT 
#community 1 
ora_1 = ora_communities_filt[ora_communities_filt['Community'] =='1'].reset_index()
#removing code R-HSA
ora_1['Term'] = ora_1['Term'].str.replace(r' R-HSA-\d+', '', regex=True)
#transform 
ora_1['neg_log10_FDR'] = -np.log10(ora_1['Adjusted P-value'])
# Overlap
ora_1[['genes_in_set',
       'genes_total']] = ora_1['Overlap'].str.split('/', expand=True).astype(int)
# Top 10 
top_1 = ora_1.nlargest(17, 'neg_log10_FDR')
top_1.drop([1,11], axis = 0, inplace = True)

plt.figure(figsize=(12,8), dpi = 300)
plt.rcParams['font.family'] = 'Arial'
sns.scatterplot(
    data=top_1,
    x='Combined Score',
    y='Term',
    size='genes_in_set',
    hue='neg_log10_FDR',
    palette='flare',
    sizes=(100,500),
    legend='brief')

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.xlabel('Combined Score',
           weight ='bold',
           labelpad = 10,
           fontsize = 14)
plt.ylabel('')
plt.title('ORA Enriched Pathways', 
          fontsize = 16,
          weight ='bold',
          pad = 10)

legend = plt.legend(
    bbox_to_anchor=(1.10, 1),
    loc='upper left',
    frameon=False,
    prop={'size': 14})

for text in legend.get_texts():
    label = text.get_text()
    if 'neg_log10_FDR' in label:
        text.set_text(r'$\mathbf{-Log_{10}(FDR)}$')  
    elif 'genes_in_set' in label:
        text.set_text(r'$\mathbf{Genes\ in\ Set}$') 

plt.tight_layout()
plt.show()

#community 14 
ora_14 = ora_communities_filt[ora_communities_filt['Community'] =='14'].reset_index()
#removing code R-HSA
ora_14['Term'] = ora_14['Term'].str.replace(r' R-HSA-\d+', '', regex=True)
#transform 
ora_14['neg_log10_FDR'] = -np.log10(ora_14['Adjusted P-value'])
# Overlap
ora_14[['genes_in_set',
       'genes_total']] = ora_14['Overlap'].str.split('/', expand=True).astype(int)
# Top 10 
top_14 = ora_14.nlargest(16, 'neg_log10_FDR')
top_14.drop(2,  axis = 0, inplace = True)

plt.figure(figsize=(12,8), dpi = 300)
plt.rcParams['font.family'] = 'Arial'
sns.scatterplot(
    data=top_14,
    x='Combined Score',
    y='Term',
    size='genes_in_set',
    hue='neg_log10_FDR',
    palette='flare',
    sizes=(100,500),
    legend='brief')

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.xlabel('Combined Score',
           weight ='bold',
           labelpad = 10,
           fontsize = 14)
plt.ylabel('')
plt.title('ORA Enriched Pathways', 
          fontsize = 16,
          weight ='bold',
          pad = 10)

legend = plt.legend(
    bbox_to_anchor=(1.10, 1),
    loc='upper left',
    frameon=False,
    prop={'size': 14})

for text in legend.get_texts():
    label = text.get_text()
    if 'neg_log10_FDR' in label:
        text.set_text(r'$\mathbf{-Log_{10}(FDR)}$')  
    elif 'genes_in_set' in label:
        text.set_text(r'$\mathbf{Genes\ in\ Set}$') 

plt.tight_layout()
plt.show()


#community 20
ora_20 = ora_communities_filt[ora_communities_filt['Community'] =='20'].reset_index()
#removing code R-HSA
ora_20['Term'] = ora_20['Term'].str.replace(r' R-HSA-\d+', '', regex=True)
#transform 
ora_20['neg_log10_FDR'] = -np.log10(ora_20['Adjusted P-value'])
# Overlap
ora_20[['genes_in_set',
       'genes_total']] = ora_20['Overlap'].str.split('/', expand=True).astype(int)
# Top 15 
top_20 = ora_20.nlargest(16, 'neg_log10_FDR')
top_20.drop(12, axis =0, inplace = True)
plt.figure(figsize=(12,8), dpi = 300)
plt.rcParams['font.family'] = 'Arial'
sns.scatterplot(
    data=top_20,
    x='Combined Score',
    y='Term',
    size='genes_in_set',
    hue='neg_log10_FDR',
    palette='flare',
    sizes=(100,500),
    legend='brief')

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.xlabel('Combined Score',
           weight ='bold',
           labelpad = 10,
           fontsize = 14)
plt.ylabel('')
plt.title('ORA Enriched Pathways', 
          fontsize = 16,
          weight ='bold',
          pad = 10)

legend = plt.legend(
    bbox_to_anchor=(1.10, 1),
    loc='upper left',
    frameon=False,
    prop={'size': 14})

for text in legend.get_texts():
    label = text.get_text()
    if 'neg_log10_FDR' in label:
        text.set_text(r'$\mathbf{-Log_{10}(FDR)}$')  
    elif 'genes_in_set' in label:
        text.set_text(r'$\mathbf{Genes\ in\ Set}$') 

plt.tight_layout()
plt.show()

#%% Graphic of Networks 

##### Degree Centrality
# Filter Communities based on interess - Metabolism, inflamation and vascular
communities = ['1', '14','20']
community_filter_sub = community_filter[community_filter['community'].isin(communities)]

# Filter edges
genes_sub = set(community_filter_sub['genes'])
network_sub = network_community[
    (network_community['preferredName_A'].isin(genes_sub)) &
    (network_community['preferredName_B'].isin(genes_sub))
]

#Network
G = nx.Graph()
for _, row in community_filter_sub.iterrows():
    G.add_node(row['genes'],
               community=row['community'],
               centrality=row['degree_centrality']) ## network metrics

for _, row in network_sub.iterrows():
    G.add_edge(row['preferredName_A'],
               row['preferredName_B'],
               weight=row['score'])

# Aesthetic
color_map = {'1': 'red', '14': 'orange', '20': 'green'}
node_colors = [color_map[G.nodes[n]['community']] for n in G.nodes()]
node_size = [40000 * G.nodes[n]['centrality'] for n in G.nodes()]
edge_widths = [d['weight'] * 2.5 for (_, _, d) in G.edges(data=True)]

#Layout
pos = nx.kamada_kawai_layout(G, weight='weight', scale=4)

# 5 TOP Hubs
top_nodes = sorted(G.nodes(data=True),
                   key=lambda x: x[1]['centrality'], reverse=True)[:10]
labels = {n: n for n, _ in top_nodes}

#plot
plt.figure(figsize=(9, 8), dpi =300)
plt.rcParams['font.family'] = 'Arial'
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.4, edge_color='gray')

# Nodes
nx.draw_networkx_nodes(G, pos,
                       node_color=node_colors,
                       node_size=node_size,
                       alpha=0.6,
                       edgecolors='black',
                       linewidths=0.8)

#Hubs labels
for n, data in top_nodes:     
    nx.draw_networkx_labels(G, pos,
                        labels=labels,
                        font_size=10,
                        font_color='black')

legend_elements = [
    mpatches.Patch(color='red', label='Immune Response'),
    mpatches.Patch(color='orange', label='Metabolism'),
    mpatches.Patch(color='green', label='Vascular and Tissue development')]

plt.legend(handles=legend_elements,
           loc='lower right',
           fontsize=12,
           frameon=False,
           title='Functional Modules',
           title_fontsize=14)

plt.title("Functional Network of Placental Genes \nin Gestational Diabetes\n(Degree Centrality)",
          fontsize=18,fontweight='bold', pad=10)
plt.axis('off')
plt.tight_layout()
plt.show()


##### Betweenness
# Filter Communities based on interess - Metabolism, inflamation and vascular
communities = ['1', '14','20']
community_filter_sub = community_filter[community_filter['community'].isin(communities)]

# Filter edges
genes_sub = set(community_filter_sub['genes'])
network_sub = network_community[
    (network_community['preferredName_A'].isin(genes_sub)) &
    (network_community['preferredName_B'].isin(genes_sub))
]

#Network
G = nx.Graph()
for _, row in community_filter_sub.iterrows():
    G.add_node(row['genes'],
               community=row['community'],
               centrality=row['betweenness']) ## network metrics

for _, row in network_sub.iterrows():
    G.add_edge(row['preferredName_A'],
               row['preferredName_B'],
               weight=row['score'])

# Aesthetic
color_map = {'1': 'red', '14': 'orange', '20': 'green'}
node_colors = [color_map[G.nodes[n]['community']] for n in G.nodes()]
node_size = [40000 * G.nodes[n]['centrality'] for n in G.nodes()]
edge_widths = [d['weight'] * 2.5 for (_, _, d) in G.edges(data=True)]

#Layout
pos = nx.kamada_kawai_layout(G, weight='weight', scale=4)

# 5 TOP Hubs
top_nodes = sorted(G.nodes(data=True),
                   key=lambda x: x[1]['centrality'], reverse=True)[:10]
labels = {n: n for n, _ in top_nodes}

#plot
plt.figure(figsize=(9, 8), dpi =300)
plt.rcParams['font.family'] = 'Arial'
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.4, edge_color='gray')

# Nodes
nx.draw_networkx_nodes(G, pos,
                       node_color=node_colors,
                       node_size=node_size,
                       alpha=0.6,
                       edgecolors='black',
                       linewidths=0.8)

#Hubs labels
for n, data in top_nodes:     
    nx.draw_networkx_labels(G, pos,
                        labels=labels,
                        font_size=10,
                        font_color='black')

legend_elements = [
    mpatches.Patch(color='red', label='Immune Response'),
    mpatches.Patch(color='orange', label='Metabolism'),
    mpatches.Patch(color='green', label='Vascular and Tissue development')]

plt.legend(handles=legend_elements,
           loc='lower right',
           fontsize=12,
           frameon=False,
           title='Functional Modules',
           title_fontsize=14)

plt.title("Functional Network of Placental Genes \nin Gestational Diabetes\n(Betweenness Centrality)",
          fontsize=18,fontweight='bold', pad=10)
plt.axis('off')
plt.tight_layout()
plt.show()

#%% Putting together the metrics 

#Communities 1 14 20
select_comm = community_filter[
    community_filter['community'].isin(['1','14','20'])]

#FoldChage in DGE
select_degs = degs[degs['Symbol'].isin(select_comm['genes'])]
select_degs = select_degs[['logFC', 'AveExpr', 'P.Value','Symbol']]

#merging metrics
all_metrics = pd.merge(select_comm, select_degs, left_on='genes',
                       right_on ='Symbol', how= 'inner')
all_metrics.drop('Symbol', axis = 1, inplace =True)

#%% Fold Change in the network

communities = ['1', '14', '20']
community_filter_sub = all_metrics[all_metrics['community'].astype(str).isin(communities)]

genes_sub = set(community_filter_sub['genes'])
network_sub = network_community[
    (network_community['preferredName_A'].isin(genes_sub)) &
    (network_community['preferredName_B'].isin(genes_sub))
]

G = nx.Graph()
for _, row in community_filter_sub.iterrows():
    G.add_node(row['genes'],
               community=str(row['community']),
               centrality=row['degree_centrality'],
               logFC=row['logFC'])

for _, row in network_sub.iterrows():
    G.add_edge(row['preferredName_A'],
               row['preferredName_B'],
               weight=row['score'])


color_map = {'1': 'red', '14': 'orange', '20': 'green'}
node_colors = [color_map[G.nodes[n]['community']] for n in G.nodes()]
node_size = [40000 * G.nodes[n]['centrality'] for n in G.nodes()]
edge_widths = [d['weight'] * 2.5 for (_, _, d) in G.edges(data=True)]


pos = nx.kamada_kawai_layout(G, weight='weight', scale=4)

top_fc_nodes = sorted(G.nodes(data=True),
                      key=lambda x: abs(x[1]['logFC']),
                      reverse=True)[:10]
labels = {n: n for n, _ in top_fc_nodes}


plt.figure(figsize=(9, 8), dpi=300)
plt.rcParams['font.family'] = 'Arial'

nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.4, edge_color='gray')

nx.draw_networkx_nodes(G, pos,
                       node_color=node_colors,
                       node_size=node_size,
                       alpha=0.6,
                       edgecolors='black',
                       linewidths=0.8)

for n, data in top_fc_nodes:
    nx.draw_networkx_labels(G, pos,
                            labels=labels,
                            font_size=12,
                            font_color='black')


legend_elements = [
    mpatches.Patch(color='red', label='Immune Response'),
    mpatches.Patch(color='orange', label='Metabolism'),
    mpatches.Patch(color='green', label='Vascular and Tissue development')
]

plt.legend(handles=legend_elements,
           loc='lower right',
           fontsize=12,
           frameon=False,
           title='Functional Modules',
           title_fontsize=14)

plt.title("Functional Network of Placental Genes \nin Gestational Diabetes\n(FoldChange)",
          fontsize=18, fontweight='bold', pad=10)
plt.axis('off')
plt.tight_layout()
plt.show()



#%% Scatter - FC x communites

#1
community_1 = all_metrics[all_metrics['community'] == '1'].sort_values('logFC', ascending = False)
degree_1 = community_1['degree_centrality'].round(3)
plt.figure(figsize=(8, 8), dpi =300)
plt.rcParams['font.family'] = 'Arial'
sns.scatterplot(
    data= community_1,
    x='logFC', 
    y='genes', 
    hue = degree_1,      
    palette='flare',     
    s=200,                   
    alpha=0.7)
plt.xlim(-2 , 2.5)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 14)
plt.ylabel("Genes", 
           fontsize = 16,
           fontweight ='bold')
plt.xlabel('Log(FoldChange)',
           fontsize = 16,
           fontweight ='bold')
plt.title('Community:Immune Response \nGDM vs Control',
          fontsize = 18,
          fontweight ='bold')
plt.legend(
           loc='lower right',
           fontsize=12,
           frameon=False,
           title='Degree Centrality',
           title_fontsize=14)

sns.despine()
plt.tight_layout()
plt.show()

#14    
community_14 = all_metrics[all_metrics['community'] == '14'].sort_values('logFC', ascending = False)
degree_14 = community_14['degree_centrality'].round(3)

plt.figure(figsize=(8, 8), dpi =300)
plt.rcParams['font.family'] = 'Arial'
sns.scatterplot(
    data= community_14,
    x='logFC', 
    y='genes', 
    hue =degree_14,      
    palette='flare',     
    s=200,                   
    alpha=0.7)
plt.xlim(-2 , 2.5)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 14)
plt.ylabel("Genes", 
           fontsize = 16,
           fontweight ='bold')
plt.xlabel('Log(FoldChange)',
           fontsize = 16,
           fontweight ='bold')
plt.title('Community: Metabolism \nGDM vs Control',
          fontsize = 18,
          fontweight ='bold')
plt.legend(
           loc='lower right',
           fontsize=12,
           frameon=False,
           title='Degree Centrality',
           title_fontsize=14)

sns.despine()
plt.tight_layout()
plt.show()

#20
community_20 = all_metrics[all_metrics['community'] == '20'].sort_values('logFC', ascending = False)
degree_20 = community_20['degree_centrality'].round(3)

plt.figure(figsize=(8, 8), dpi =300)
plt.rcParams['font.family'] = 'Arial'
sns.scatterplot(
    data= community_20,
    x='logFC', 
    y='genes', 
    hue = degree_20,      
    palette='flare',     
    s=200,                   
    alpha=0.7)
plt.xlim(-2 , 2.5)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 14)
plt.ylabel("Genes", 
           fontsize = 16,
           fontweight ='bold')
plt.xlabel('Log(FoldChange)',
           fontsize = 16,
           fontweight ='bold')
plt.title('Community: Vascular and Tissue development \nGDM vs Control',
          fontsize = 18,
          fontweight ='bold')
plt.legend(
           loc='lower right',
           fontsize=12,
           frameon=False,
           title='Degree Centrality',
           title_fontsize=14)

sns.despine()
plt.tight_layout()
plt.show()



