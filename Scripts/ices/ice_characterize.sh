#!/bin/bash

# Default paths
DB_PATH="/media/usuario/datos/Databases"
PLASMIDFINDER_DB="$DB_PATH/orit/oriT.fasta"
INTEGRASES_HMM="$DB_PATH/integrases_hmm/integrases_ser_try.hmm"
PLASMIDFINDER_REF="$DB_PATH/plasmidfinder/plasmidfinder.fna"
BAKTA_DB="$DB_PATH/bakta/db"

# Default number of threads
THREADS=30

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_file> [-o <output_dir>] [-t <threads>] [-h]"
}

# Parse command line options
while getopts "i:o:t:h" opt; do
    case $opt in
        i)
            INPUT_FILE="$OPTARG"
            ;;
        o)
            OUTPUT_DIR="$OPTARG"
            ;;
        t)
            THREADS="$OPTARG"
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            usage
            exit 1
            ;;
    esac
done

# Check for required input file
if [ -z "$INPUT_FILE" ]; then
    echo "Error: Input file is required."
    usage
    exit 1
fi

# Set default output directory to the directory of the input file
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$(dirname "$INPUT_FILE")/anotacion
fi

# Get the directory of the input file
INPUT_DIR=$(dirname "$INPUT_FILE")

# Create output directories
mkdir -p "$OUTPUT_DIR"/bakta
mkdir -p "$OUTPUT_DIR"/conj
mkdir -p "$OUTPUT_DIR"/txss
mkdir -p "$OUTPUT_DIR"/MOBscan
mkdir -p "$OUTPUT_DIR"/oriT
mkdir -p "$OUTPUT_DIR"/rep
mkdir -p "$OUTPUT_DIR"/integrones
mkdir -p "$OUTPUT_DIR"/integrasas

# Load conda function
source /home/usuario/miniconda3/etc/profile.d/conda.sh

# Iterate over each line in the input file
# while IFS= read -r i; do
while read -r i; do
    f_name=$(basename "$INPUT_DIR"/"$i" .fasta)
    echo $f_name
    # Activate icecreen environment
    # conda activate nano2

    # Annotate with Bakta
    bakta --db "$BAKTA_DB" --keep-contig-headers --skip-plot --verbose --output "$OUTPUT_DIR"/bakta --threads "$THREADS" --force "$INPUT_DIR"/"$i"
    
    # Run macsyfinder for CONJScan
    mkdir -p "$OUTPUT_DIR"/conj/"$f_name"
    macsyfinder --force -w $THREADS --timeout 30m -m CONJScan/Plasmids all --sequence-db "$OUTPUT_DIR"/bakta/"$f_name".faa --db-type ordered_replicon -o "$OUTPUT_DIR"/conj/"$f_name" -vv --e-value-search 0.0001 

    # Run macsyfinder for TXSScan
    mkdir -p "$OUTPUT_DIR"/txss/"$f_name"
    macsyfinder --force -w $THREADS --timeout 30m -m TXSScan all --sequence-db "$OUTPUT_DIR"/bakta/"$f_name".faa --db-type ordered_replicon -o "$OUTPUT_DIR"/txss/"$f_name" -vv --e-value-search 0.0001
    
    # Run MOBscan
    mkdir -p "$OUTPUT_DIR"/MOBscan
    mobscan.py -n 28 "$OUTPUT_DIR"/bakta/"$f_name".faa "$OUTPUT_DIR"/MOBscan/"$f_name"
    
    # Run blastn for oriT
    mkdir -p "$OUTPUT_DIR"/oriT
    blastn -query "$INPUT_DIR"/"$i" -db "$PLASMIDFINDER_DB" -outfmt "6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send qframe evalue bitscore" -perc_identity 90 -max_target_seqs 1 | sort -rk2 | uniq > "$OUTPUT_DIR"/oriT/"$f_name".csv
    
    # Run blastn for rep Plasmid
    mkdir -p "$OUTPUT_DIR"/rep
    blastn -query "$INPUT_DIR"/"$i" -db "$PLASMIDFINDER_REF" -outfmt "6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send qframe evalue bitscore" -perc_identity 90 -max_target_seqs 1 | sort -rk2 | uniq > "$OUTPUT_DIR"/rep/"$f_name".csv
    
    # Activate integron_finder environment
    conda activate integron_finder
    
    # Run integron_finder
    integron_finder "$INPUT_DIR"/"$i" --cpu "$THREADS" --outdir "$OUTPUT_DIR"/integrones
    
    # Deactivate integron_finder environment
    conda deactivate
    
    # Run hmmsearch for integrases
    mkdir -p "$OUTPUT_DIR"/integrasas
    hmmsearch --tblout "$OUTPUT_DIR"/integrasas/"$f_name".tab "$INTEGRASES_HMM" "$OUTPUT_DIR"/bakta/"$f_name".faa
    sed -i -E 's/\s{2,}/\t/g' "$OUTPUT_DIR"/integrasas/"$f_name".tab

    # Run abricate for AbR and virulence genes
    mkdir -p "$OUTPUT_DIR"/abr
    abricate --minid 75 "$INPUT_DIR"/"$i" > "$OUTPUT_DIR"/abr/"$f_name".tab
    abricate --summary "$OUTPUT_DIR"/abr/"$f_name".tab > "$OUTPUT_DIR"/abr/"$f_name".summary

done < "$INPUT_FILE"
