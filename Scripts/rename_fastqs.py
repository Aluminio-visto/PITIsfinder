# %%
import pandas as pd
import os

# %% [markdown]
# Datos a cambiar

# %%
# Este archivo tiene que tener la fecha en formato YYYY_MM_DD y debe coincidir con el nombre de la carpeta donde se almacenan los fastqs
infile = '/media/usuario/datos2/Mepram/repositorio2/datos_seq_nuevos.csv'
df = pd.read_csv(infile, dtype={'Barcode': str, 'BC (rep)': str, 'BC (rep2)': str})

# Directorio madre
inpath = '/media/usuario/datos2/Mepram/'
# Directorio output
outpath = '/media/usuario/datos2/Mepram/repositorio2/01_reads/'

# %% [markdown]
# Ejecución

# %%
# Filtramos columna sin ID único
df.dropna(subset=['ID único'], inplace=True)

# Filtramos columna con ID incompleto
df = df.loc[~df['ID único'].str.contains('__')]
df = df.loc[~df['ID único'].str.contains(r'\?')]

# fill NAs
df.fillna('', inplace=True)

# %%
for _, row in df.iterrows():
    id_unico = row['ID único']
    bc1 = f"{row['Barcode']}"
    bc2 = f"{row['BC (rep)']}"
    bc3 = f"{row['BC (rep2)']}"
    date1 = row['Fecha seq']
    date2 = row['Fecha seq (rep)']
    date3 = row['Fecha seq (rep2)']
    # print(id_unico)

    if(bc3):
        print(id_unico, f'{inpath}{date1}/fastq_pass/barcode{int(bc1):02d}/*', f'{inpath}{date2}/fastq_pass/barcode{int(bc2):02d}/*', f'{inpath}{date3}/fastq_pass/barcode{int(bc3):02d}/*')
        cat_cmd = ['cat', f'{inpath}{date1}/fastq_pass/barcode{int(bc1):02d}/*', f'{inpath}{date2}/fastq_pass/barcode{int(bc2):02d}/*', f'{inpath}{date3}/fastq_pass/barcode{int(bc3):02d}/*', '>', f'{outpath}/{id_unico}.fastq.gz']
        os.system(" ".join(cat_cmd))
    elif(bc2):
        print(id_unico, f'{inpath}{date1}/fastq_pass/barcode{int(bc1):02d}/*', f'{inpath}{date2}/fastq_pass/barcode{int(bc2):02d}/*')
        cat_cmd = ['cat', f'{inpath}{date1}/fastq_pass/barcode{int(bc1):02d}/*', f'{inpath}{date2}/fastq_pass/barcode{int(bc2):02d}/*', '>', f'{outpath}/{id_unico}.fastq.gz']
        os.system(" ".join(cat_cmd))
    elif(date1):
        print(id_unico, f'{inpath}{date1}/fastq_pass/barcode{int(bc1):02d}/*')
        cat_cmd = ['cat', f'{inpath}{date1}/fastq_pass/barcode{int(bc1):02d}/*', '>', f'{outpath}/{id_unico}.fastq.gz']
        os.system(" ".join(cat_cmd))
