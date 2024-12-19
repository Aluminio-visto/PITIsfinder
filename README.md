# PITIsfinder
## Bacterial whole-genome analysis pipeline for Mobile Genetic Element Detection and Classification

This pipeline is currently in construction. Translating to Next-Flow the whole pipeline. 
You may find some input-output examples in the folder examples.
Further installation instructions will be added when environments are solved. 
You may take a look to the main scripts in Scripts folder and main.sh.
The virtual conda environments used are in the environemnts folder.

![image](https://github.com/user-attachments/assets/b36a0998-3af2-44d6-b167-e54369af821c)

## Input
The working directory must contain the following files and folders/subfolders:
**Sample input, output and intermediate files can be found at: https://zenodo.org/records/14523202**
- fastq_pass/barcode*/*.fastq.gz
- final_summary_*.txt
- report_*.json
- lista_seq.tsv
- datos_analisis.csv
- datos_seq.csv

## Output
**Sample input, output and intermediate files can be found at: https://zenodo.org/records/14523202**
```
.
├── 01_reads/
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   ├── ...
│   ├── QC/
│   │   ├── sample1/
│   │   │   └── html, png, txt and log files from Nanoplot
│   │   └── ...
│   ├── QC_reads.csv
│   └── QC_reads_pre.csv
├── 02_filter/
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   ├── ...
│   ├── QC/
│   │   ├── sample1/
│   │   │   └── html, png, txt and log files from Nanoplot
│   │   └── ...
│   ├── QC_reads.csv
│   └── QC_reads_post.csv
├── 03_assemblies/
│   ├── sample1/
│   │   ├── assembly.fasta
│   │   ├── assembly_graph.gfa
│   │   ├── deconcat/
│   │   │   ├── bam, png and fasta files
│   │   │   └── ...
│   │   └── ...
│   ├── sample.../
│   └── quast/
├── 04_taxonomies/
│   ├── kraken2/
│   │   ├── sample1.out
│   │   ├── sample1.report
│   │   ├── ...
│   │   ├── genus.csv
│   │   ├── species.csv
│   │   └── strains.csv
│   ├── ectyper/
│   │   ├── blast_output_alleles.txt
│   │   └── output.tsv
│   ├── gambit.csv
│   ├── gtdb/
│   |   ├── sample1.fasta
│   |   ├── sample2.fasta
│   |   └── ...
│   └── kraken_report.csv
├── 05_plasmids/
│   ├── plasmid1_sample1.fasta
│   ├── plasmid2_sample1.fasta
│   └── ...
├── 08_Anotacion/
│   ├── sample1/
│   │   ├── abricate/
│   │   │   └── sample1.tab
│   │   ├── copla/
│   │   │   ├── fasta and tsv files
│   │   ├── mob_recon/
│   │   │   ├── replicons fasta files
│   │   │   ├── contig_report.txt
│   │   │   └── mobtyper_results.txt
│   │   └── Bakta annotation files (gff, gbk, faa, tsv...)
│   ├── sample.../
│   └── AbR.tab
├── 09_phages/
│   ├── phastest_deep/
│   │   ├── sample1/
│   │   │   ├── PHASTEST output files
│   │   └── sample.../
│   ├── PHAGE_phageID_sample1.fasta
│   └── phage_summary.csv
├── 11_integrons/
│   ├── sample1/
│   │   └── Results_Integron_Finder_Eclo_VC_28100/
│   │   │   ├── Integron Finder output files
│   │   │   ├── Abricate output files
│   │   │   └── Prokka output files
│   ├── sample.../
│   ├── int_cassette1_cassette2_sample1_1.fasta
│   └── integron_summary.csv
├── QC_reads.csv
├── QC_assembly.csv
├── AbR_modif.csv
├── taxonomy.csv
├── mlst_modif.csv
├── kleborate.tsv
├── copla_modif.csv
├── datos_analisis_nuevos.csv
└── datos_seq_nuevos.csv                     
```

## Execution
Before the initial execution, make sure all conda environments available at the environments folder are installed and all Databases needed in the pipeline are locally accessible in the specified locations at main.sh.

Sequentially execute the commands specified at main.sh.
