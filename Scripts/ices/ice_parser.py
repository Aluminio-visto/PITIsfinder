#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
import re
import math
import csv
import time
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import matplotlib.pyplot as plt
import argparse

# # Parse output

# Caracterizados:
# 
# - tRNA (Bakta/tRNAscan-SE)
# - Integrasas (Bakta, integron_finder, HMM de Ser/Try recombinasas)
# - oriT (Blastn)
# - Relaxosoma (macsyfinder CONJScan, MOBscan)
# - T4SS (macsyfinder TXSS).
# - Rep plásmidos (Blastn)
# - att (integron_finder)
# - Cargo (Bakta, AbR?)    

# In[2]:

parser = argparse.ArgumentParser(description='ICE Parser Script',
                                 epilog='Example: python ice_parser.py -p /path/to/your/data')
parser.add_argument('-p', '--path', type=str, help='Path to data directory')
parser.add_argument('-n', '--ice_nts', type=int, default=30000, help='Upstream and downstream ICE length (default: 30000)')


# Parse command-line arguments
args = parser.parse_args()
path = args.path
ICE_NTS = args.ice_nts

# Load sample names
with open(f'{path}samples.txt') as in_f:
    samples = [".".join(fasta.split(".")[:-1]) for fasta in in_f.read().splitlines()]


# In[3]:


def extract_bakta(bakta_path):
    # Bakta (anota con tRNAscan-SE)
    
    # Read bakta annotation file
    if os.path.exists(bakta_path) and os.path.getsize(bakta_path) > 0:
        df_bakta = pd.read_csv(bakta_path, sep='\t', comment='#', header=None)
        df_bakta.columns = ['Sequence Id', 'Type', 'Start', 'Stop', 'Strand', 'Locus Tag', 'Gene', 'Product', 'DbXrefs']
        df_bakta['Type'].replace('cds', 'CDS', inplace=True)
        df_bakta['Tag'] = 'Cargo'
        try:
            # Keep relevant features(tRNA, integrase...) <- DE MOMENTO ESTO, SI SE NOS OCURRE ALGO MÁS, BIENVENIDO SEA
            df_trna = df_bakta[df_bakta['Type'] == 'tRNA']
            df_trna = df_trna[['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product']]
            # df_trna = df_trna[df_trna['Product'].str.lower().str.contains('trna-ser|trna-tyr')]
            df_trna['Source'] = 'Bakta'
            df_trna['Sys_id'] = 'tRNA'
            df_trna['Tag'] = 'tRNA'
        except:
            df_trna = None

        try:
            df_integrase = df_bakta.copy()
            df_integrase = df_integrase[df_integrase['Product'].str.lower().str.contains('integrase')]
            df_integrase = df_integrase[['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product']]
            df_integrase['Sys_id'] = 'Integrase'
            df_integrase['Source'] = 'Bakta'
            df_integrase['Tag'] = 'Integrase'
        except:
            df_integrase = None
     
        return df_bakta, df_trna, df_integrase

    else: return None, None, None


# In[4]:


def extract_orit(orit_path):

    if os.path.exists(orit_path) and os.path.getsize(orit_path) > 0:
        df_orit = pd.read_csv(orit_path, sep='\t' , header=None)
        df_orit.columns = ['Sequence Id', 'Gene', 'pident', 'qcovhsp', 'length', 'qlen', 'slen', 'Start', 'Stop', 'sstart', 'send', 'qframe', 'evalue', 'bitscore']
        ######QUE EN UN FUTURO DIGA LA STRAND???######
        df_orit['Strand'] = pd.NA
        df_orit['Product'] = df_orit['Gene']
        df_orit = df_orit[['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product']]
        df_orit['Sys_id'] = 'OriT'
        df_orit['Source'] = 'Blastn'
        df_orit['Tag'] = 'OriT'
        return df_orit

    else: return None


# In[5]:


def extract_rep(rep_path):

    if os.path.exists(rep_path) and os.path.getsize(rep_path) > 0:
        df_rep = pd.read_csv(rep_path, sep='\t' , header=None)
        df_rep.columns = ['Sequence Id', 'Gene', 'pident', 'qcovhsp', 'length', 'qlen', 'slen', 'Start', 'Stop', 'sstart', 'send', 'qframe', 'evalue', 'bitscore']
        ######QUE EN UN FUTURO DIGA LA STRAND???######
        df_rep['Strand'] = pd.NA
        df_rep['Product'] = df_rep['Gene']
        df_rep = df_rep[['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product']]
        df_rep['Sys_id'] = 'Rep'
        df_rep['Source'] = 'Blastn'
        df_rep['Tag'] = 'Rep'
        return df_rep

    else: return None


# In[6]:


def extract_mob(mob_path, df_bakta):
    if os.path.exists(mob_path) and not df_bakta.empty:
        try:
            df_mob = pd.read_csv(mob_path, sep='\t', comment='#')
            df_mob = pd.merge(df_bakta, df_mob, left_on='Locus Tag', right_on='hit_id')
            df_mob[['Product', 'Gene']] = df_mob['gene_name'].str.extract('(.*[_\w]?)_(\w+)')
            df_mob['Sys_id'] = df_mob['sys_id'].str.replace(f'{s}_', '')
            df_mob = df_mob[['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product', 'Sys_id']]
            
            # Create a new column 'Sys_id_sub' with updated values, divides extra large (presumably wrong) 
            df_mob['Sys_id_sub'] = df_mob.groupby((df_mob['Start'] - 1) // 50000).grouper.group_info[0] + 1
            df_mob['Sys_id'] = df_mob['Sys_id_sub'].apply(lambda x: f"{df_mob['Sys_id'].iloc[0]}_{x}")

            df_mob['Source'] = 'CONJScan'
            df_mob['Tag'] = 'T4SS'
            df_mob.drop(columns='Sys_id_sub', inplace=True)
            return df_mob
        except: pass
    return None


# In[7]:


def extract_mobscan(mobscan_path, df_bakta):
    if os.path.exists(mobscan_path) and not df_bakta.empty:
        try:
            df_mobscan = pd.read_csv(mobscan_path, sep='\t', comment='#')
            df_mobscan = pd.merge(df_bakta, df_mobscan, left_on='Locus Tag', right_on='Query name')
            df_mobscan[['Product_y', 'Gene_y']] = df_mobscan['Profile HMM'].str.extract('(\w+)_(\w+)')
            df_mobscan = df_mobscan[['Sequence Id', 'Start_x', 'Stop', 'Strand', 'Product_y', 'Gene_y']]
            df_mobscan = df_mobscan.rename(columns={'Start_x': 'Start',
            'Product_y': 'Product',
                                                    'Gene_y': 'Gene'})
            df_mobscan['Source'] = 'MOBscan'
            df_mobscan['Tag'] = 'MOB'
            return df_mobscan
        except: pass
    return None


# In[8]:


def extract_conj(conj_path, df_bakta):
    if os.path.exists(conj_path) and not df_bakta.empty:
        try:
            df_conj = pd.read_csv(conj_path, sep='\t', comment='#')
            df_conj = pd.merge(df_bakta, df_conj, left_on='Locus Tag', right_on='hit_id')
            df_conj[['Product', 'Gene']] = df_conj['gene_name'].str.extract('(.*[_\w]?)_(\w+)')
            df_conj = df_conj[df_conj['Product'] == 'T4SS']
            df_conj['Sys_id'] = df_conj['sys_id'].str.replace(f'{s}_', '')
            df_conj = df_conj[['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product', 'Sys_id']]
            df_conj['Source'] = 'TXSScan'
            df_conj['Tag'] = 'T4SS'
            return df_conj
        except: pass
    return None


# In[9]:


def extract_integron(integron_path):
    try:
        df_int = pd.read_csv(integron_path, sep='\t', comment='#')
        df_int = df_int[(df_int['annotation'] == 'attC') | (df_int['annotation'] == 'intI')]
        df_int = df_int[['ID_replicon', 'pos_beg', 'pos_end', 'strand', 'annotation', 'model']]
        df_int.columns = ['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product']
        df_int['Strand'] = np.where(df_int['Strand'] == -1, '-', '+')
        df_int['Sys_id'] = df_int['Gene']
        df_int['Sys_id'].replace({'intI': 'Integrase'}, inplace = True)
        df_int['Source'] = 'Integron_finder'
        df_int['Tag'] = ['Integrase' if x == 'intI' else 'att' for x in df_int['annotation']]
        return df_int
    except:
        return None


# In[10]:


def extract_integrase(integron_path, df_bakta):
    try:
        df_int = pd.read_csv(integrase_path, sep='\t', comment='#', header=None)
        df_int.columns = ['target name', 'accession_t', 'query name', 'accession_q', 'E-value_fs', 'score_fs', 'bias_fs', 'E-value_bd', 'score_bd', 'bias_bd', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'other']
        df_int[['inc', 'description of target']] = df_int['other'].str.split(' ', n=1, expand=True)
        df_int.drop(columns=['other'], inplace=True)
        df_int = pd.merge(df_bakta, df_int, left_on='Locus Tag', right_on='target name')
        df_int = df_int.groupby('Locus Tag', group_keys=False).apply(lambda x: x.loc[x.score_fs.idxmax()])
        df_int = df_int[['Sequence Id', 'Start', 'Stop', 'Strand', 'query name', 'Product']]
        df_int['Sys_id'] = 'Integrase'
        df_int['Source'] = 'HMM_integrase'
        df_int.rename(columns={'query name': 'Gene'}, inplace=True)
        df_int['Tag'] = 'Integrase'
        return df_int
    except:
        return None



def extract_abr(abr_path):
    try:
        df_abr = pd.read_csv(abr_path, sep='\t')
        df_abr.columns = ['file','Sequence Id', 'Start', 'Stop', 'Gene', 'cov', 'cov_map', 'gaps', 'perc_cov', 'perc_id', 'db', 'acc']
        df_abr.drop(df_abr.loc[df_abr['perc_cov'] < 75].index, inplace=True)
        df_abr.drop(columns=['file', ], inplace=True)
        df_abr['Strand'] = '+'
        df_abr['Product'] = df_abr['Gene']
        df_abr = df_abr[['Sequence Id', 'Start', 'Stop', 'Strand', 'Gene', 'Product']]
        df_abr['Sys_id'] = 'AbR'
        df_abr['Source'] = 'abricate'
        df_abr['Tag'] = 'AbR'
        return df_abr
    except:
        return None

# In[11]:


def extract_type(mobrecon_path, complete_df):
    if os.path.exists(mobrecon_path) and not complete_df.empty:
        try:
            df_mobrecon = pd.read_csv(mobrecon_path, sep='\t')
            complete_df = complete_df.merge(df_mobrecon[['molecule_type', 'contig_id']], 
                        left_on='Sequence Id', right_on='contig_id', 
                                            suffixes=(False, False)).drop(columns = ['contig_id'])
            complete_df_type = complete_df[['Sequence Id', 'molecule_type', 'Start', 'Stop', 'Strand', 'Gene', 'Product', 'Sys_id', 'Source']]
            return complete_df_type, df_mobrecon
        except: pass
    return None, None


# In[12]:


def merge_rows(series):
    unique_values = [i for i in series.unique() if not(pd.isnull(i)) == True]
    if len(unique_values) == 0:
        return None
    elif len(unique_values) == 1:
        return unique_values[0]
    else:
        return ", ".join(unique_values)


# In[13]:


def get_all_sys(column):
    sys = set()
    for i in column.dropna():
        elems = (i.split(', '))
        sys.update(elems)
    filtered_sys = [i for i in sys if 'CONJ' in i or 'T4SS' in i or 'MOB' in i]
    return filtered_sys


# In[14]:


def extract_systems(systems, d_df):
    d_sys = {}
    for system in systems:
        system_df = d_df.loc[d_df['Sys_id'].str.contains(system, na=False)]
        start = system_df.iloc[0]['Start']
        end = system_df.iloc[-1]['Stop']
        d_sys[system] = (start, end)
    try:
        d_sys_sort = sorted(d_sys.items(), key=lambda x: x[1][0])
        return d_sys_sort
    except: return d_sys


# In[15]:


def extract_fasta(path, s, ice_n, k, start_sys, end_sys):
    # fasta_file = f'{path}/{s}.fasta'
    fasta_file = s + '.fasta'
    # outfile = f'{path}/fastas/ICE_{ice_n}_{k}_{s}.fasta'
    outfile = f'{path}/fastas/ICE_{ice_n}_{k}_{os.path.splitext(os.path.basename(s))[0]}.fasta'
    if start_sys < 0:
        start_sys = 0

    with open(outfile, 'w') as f_out:
        with open(fasta_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if k in record.id:
                    ice_id = f'ICE_{ice_n}_{k}_{s}'
                    if len(record.seq) >= end_sys+1:
                        ice_fasta = f'{record.seq[start_sys:end_sys+1]}'
                    else:
                        ice_fasta = f'{record.seq[start_sys:]}'
                        end_sys = len(record.seq)
                    f_out.write(f'>{ice_id}\n')
                    f_out.write(ice_fasta)
                    return ice_id, ice_fasta, start_sys, end_sys
    return None, None, None, None


# In[16]:


def characterize_ice(ice_df, ice_n, k, s, length_sys, start_sys, end_sys, sys_i_id, sys_j_id, path, df_bakta, contig_df):
    # T4CP
    try:
        t4cp = ice_df['Gene'].str.lower().str.contains('t4cp').any()
    except:
        t4cp = False

    # Virb4
    try:
        virb4 = ice_df['Gene'].str.lower().str.contains('virb4').any()
    except:
        virb4 = False
    
    # MPF and MOB types
    if 'T4SS' in sys_i_id:
        mpf_type = sys_i_id.split('_')[1]
        mob_df = ice_df.loc[ice_df['Sys_id'].str.contains(sys_j_id)]
        try:
            mob_type = mob_df.loc[mob_df['Gene'].str.lower().str.contains('mob')]['Gene'].values[0]
        except:
            mob_type = None
    else:
        mpf_type = sys_j_id.split('_')[1]
        mob_df = ice_df.loc[ice_df['Sys_id'].str.contains(sys_i_id)]
        try:
            mob_type = mob_df.loc[mob_df['Gene'].str.lower().str.contains('mob')]['Gene'].values[0]
        except:
            mob_type = None

    # Cargo
    sub_df = df_bakta.loc[df_bakta['Sequence Id'] == k]
    df_cargo = sub_df.loc[(sub_df['Start'] >= start_sys - ICE_NTS) & (sub_df['Stop'] <= end_sys + ICE_NTS)]
    cargo = list(df_cargo['Product'])

    # Extract ICE fasta
    ice_id, ice_fasta, start_sys, end_sys = extract_fasta(path, sample, ice_n, k, start_sys - ICE_NTS, end_sys + ICE_NTS)

    df = pd.merge(df_bakta, contig_df, on=['Sequence Id', 'Start', 'Stop'], how='outer', suffixes=("","_cargo"))
    df = df.loc[(df['Start'] >= start_sys) & (df['Stop'] <= end_sys)]

    # Integrase
    integrase = set(list(df['Gene'].loc[df['Sys_id'].str.contains('Integrase', na=False)]))

    # Resistance
    resist = list(df['Gene'].loc[df['Sys_id'].str.contains('AbR', na=False)])

    d_ice = {'Nombre muestra': s,
             'N ICE': ice_n,
             'contig': k,
             'longitud': length_sys,
             'start': start_sys,
             'stop': end_sys,
             'MOBtype': mob_type,
             'MPFtype': mpf_type,
             'Integrase': integrase,
             'T4CP': t4cp,
             'Virb4': virb4,
             'Cargo': cargo,
             'AbR': resist}

    # Get the list of cargo columns
    cargo_columns = [col for col in df.columns if col.endswith('_cargo')]

    # Iterate through each row and copy values from cargo columns to their corresponding columns without the suffix
    for index, row in df.iterrows():
        for cargo_col in cargo_columns:
            base_col = cargo_col.replace('_cargo', '')
            if pd.notna(row[cargo_col]):
                df.at[index, base_col] = row[cargo_col]

    # Drop the cargo columns after copying the values
    df.drop(columns=cargo_columns, inplace=True)
    dtype_mapping = {"Sequence Id": str,
                     "Type": str,
                     "Start": int,
                     "Stop": int,
                     "Strand": str,
                     "Locus Tag": str,
                     "Gene": str,
                     "Product": str,
                     "DbXrefs": str,
                     "Tag": str,
                     "Source": str,
                     "Sys_id": str}
                     
    df = df.astype(dtype_mapping)
    df.to_csv(f'{path}/tables/ICE_{ice_n}_{k}_{s}.csv')
    
    # Extract ICE gbk
    ice_object = Seq(ice_fasta)
    record = SeqRecord(ice_object,
                       id = f'{ice_id}',
                       name=f'{ice_id}',
                       annotations={"molecule_type": "DNA"})
    

    d_fasta = SeqIO.to_dict(SeqIO.parse(f'{path}/anotacion/bakta/{s}.faa', "fasta"))
    for _, row in df.iterrows():
        gene_start = int(row['Start'] - (start_sys))
        gene_stop = abs(int(row['Stop'] - (start_sys)))
        if row['Strand'] == '+':
            strand = 1
        else:
            strand = -1
        try:
            feature = SeqFeature(FeatureLocation(start=gene_start, end=gene_stop),
                                type = 'gene',
                                strand = strand,
                                id = row['Gene'],
                                qualifiers = {'locus_tag': row['Locus Tag'],
                                            'db_xref': row['DbXrefs']})
            record.features.append(feature)

            feature = SeqFeature(FeatureLocation(start=gene_start, end=gene_stop),
                                type = row['Type'],
                                strand = strand,
                                id = row['Gene'],
                                qualifiers = {'locus_tag': row['Locus Tag'],
                                            'gene': row['Gene'],
                                            # 'strand': strand,
                                            'codon_start': 1,
                                            'transl_table': 11,
                                            'product': row['Product'],
                                            'db_xref': row['DbXrefs'],
                                            'tag': row['Tag'],
                                            'translation': d_fasta[row['Locus Tag']].seq})
            record.features.append(feature)
            

        
        except:
            feature = SeqFeature(FeatureLocation(start=gene_start, end=gene_stop),
                                type = row['Type'],
                                strand = strand,
                                id = row['Gene'],
                                qualifiers = {'gene': row['Gene'],
                                            # 'strand': strand,
                                            'codon_start': 1,
                                            'transl_table': 11,
                                            'product': row['Product'],
                                            'tag': row['Tag'],
                                            'db_xref': row['DbXrefs']})
            record.features.append(feature)

    output_file = open(f'{path}/gbks/ICE_{ice_n}_{k}_{s}.gbk', 'w')
    SeqIO.write(record, output_file, 'genbank')

    ice_n += 1

    return d_ice, ice_n


# In[17]:


def system_selection(sys_pairs, df_bakta, contig_df, complete_ice_df, ice_n, k, s, path):

    new_row = None
    for sys in sys_pairs:
        sys_i_id, sys_i_pos = sys[0]
        sys_j_id, sys_j_pos = sys[1]

        # Buscar parejas de T4SS (MPF type) y MOB/CONJ (MOB type)
        if ('T4SS' in sys_i_id or 'T4SS' in sys_j_id) and ('T4SS' not in sys_i_id or not 'T4SS' in sys_j_id):
            start_sys = min([sys_i_pos[0], sys_j_pos[0]])
            end_sys = max([sys_i_pos[1], sys_j_pos[1]])
            length_sys = end_sys - start_sys
            # Distancia menor de 50000pdb
            if length_sys < 50000: #and length_sys > 20000:
                # ICE df (containing MPF and MOB)
                ice_df = contig_df.loc[(contig_df['Sys_id'].str.contains(sys_i_id)) | (contig_df['Sys_id'].str.contains(sys_j_id))]

                # Tabla con características del ICE
                d_ice, ice_n = characterize_ice(ice_df, ice_n, k, s, length_sys, start_sys, end_sys, sys_i_id, sys_j_id, path, df_bakta, contig_df)
                new_row = pd.DataFrame.from_dict(d_ice, orient='index').T
                complete_ice_df = pd.concat([complete_ice_df, new_row])
    return complete_ice_df, ice_n


# In[18]:


os.system(f'mkdir -p {path}/tables')
os.system(f'mkdir -p {path}/fastas')
os.system(f'mkdir -p {path}/gbks')


# In[19]:


# Extract table per sample

complete_ice_df = pd.DataFrame(columns=['Nombre muestra', 'N ICE', 'contig', 'longitud', 'start', 'stop', 'MOBtype', 'MPFtype', 'T4CP', 'Virb4', 'Integrase', 'Cargo', 'AbR'])

for sample in samples:
    s = os.path.splitext(os.path.basename(sample))[0]
    print(sample)

    # Gral annotation
    bakta_path = f'{path}/anotacion/bakta/{s}.tsv'
    df_bakta, df_trna, df_integrase = extract_bakta(bakta_path)

    # OriT
    orit_path = f'{path}/anotacion/oriT/{s}.csv'
    df_orit = extract_orit(orit_path)

    rep_path = f'{path}/anotacion/rep/{s}.csv'
    df_rep = extract_rep(rep_path)

    # Mobilizable module
    mob_path = f'{path}/anotacion/conj/{s}/best_solution.tsv'
    df_mob = extract_mob(mob_path, df_bakta)

    # MOBscan
    mobscan_path = f'{path}/anotacion/MOBscan/{s}/results_tab.tsv'
    df_mobscan = extract_mobscan(mobscan_path, df_bakta)

    # Conjugative module
    conj_path = f'{path}/anotacion/txss/{s}/best_solution.tsv'
    df_conj = extract_conj(conj_path, df_bakta)

    # Integrons
    integron_path = f'{path}/anotacion/integrones/Results_Integron_Finder_{s}/{s}.integrons'
    df_integron = extract_integron(integron_path)

    # Integrases
    integrase_path = f'{path}/anotacion/integrasas/{s}.tab'
    df_hmm = extract_integrase(integrase_path, df_bakta)

    abr_path = f'{path}/anotacion/abr/{s}.tab'
    df_abr = extract_abr(abr_path)

    # Merge all dataframes
    complete_df = pd.concat([df_trna, df_integrase, df_orit, df_rep, df_mob, df_mobscan, df_conj, df_integron, df_hmm, df_abr], axis=0)
    complete_df.sort_values(['Sequence Id', 'Start'], inplace=True)

    # # Assign Chr/plasmid
    # mobrecon_path = f'{path}/anotacion/plasmids/{s}/contig_report.txt'
    # complete_df_type, df_mobrecon = extract_type(mobrecon_path, complete_df)

    # Merge duplicate entries called by different tools
    complete_df_nodup = complete_df.groupby(['Sequence Id', 'Start', 'Stop'], as_index=False).agg(merge_rows)
    complete_df_nodup.to_csv(f'{path}tables/{s}.csv', index=False)

    # Group by contig
    d_df = {i: x for i, x in complete_df_nodup.groupby('Sequence Id', as_index=False)}
    
    # Extract all systems in contig (contig by contig) and process
    for k in d_df.keys():
        systems = get_all_sys(d_df[k]['Sys_id'].dropna())
        # Extract systems
        d_sys = extract_systems(systems, d_df[k])
        ice_n = 1

        # Get system pairs
        sys_pairs = []
        for i in range(len(d_sys)):
            for j in range(i + 1, len(d_sys)):
                sys_pairs.append([d_sys[i], d_sys[j]])

        # Loop through systems extracting characteristics
        complete_ice_df, ice_n = system_selection(sys_pairs, df_bakta, d_df[k], complete_ice_df, ice_n, k, s, path)


# In[20]:

#['Nombre muestra', 'N ICE', 'contig', 'longitud', 'start', 'stop', 'MOBtype', 'MPFtype', 'T4CP', 'Virb4', 'Cargo', 'AbR', 'Integrase']
complete_ice_df = complete_ice_df[['Nombre muestra', 'contig', 'N ICE', 'longitud', 'start', 'stop', 'MOBtype', 'MPFtype', 'Integrase', 'AbR', 'T4CP', 'Virb4']]
complete_ice_df['Mobility'] = complete_ice_df['T4CP'] & complete_ice_df['Virb4']
complete_ice_df.to_csv(f'{path}/tables/ICE_summary.csv', index=False)
