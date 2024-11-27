rename_fastqs() {
    mkdir 01_reads
    echo ID$'\t'R1$'\t'R2$'\t'LongFastQ$'\t'Fast5$'\t'GenomeSize$'\t'Fasta > samplesheet.tsv
    tail -n +2 lista_seq.tsv | while read -r cult cep id bc c rep; do
    echo $id $bc $rep
    # Si no es repetido
    if [ -z "$rep" ]; then
        # Renombrar
        cat fastq_pass/$bc/*.fastq.gz > 01_reads/"$id".fastq.gz;
        # Actualizar repo
        cp 01_reads/"$id".fastq.gz /home/usuario/Seqs/Mepram/repositorio/01_reads/$id".fastq.gz";
    # Si es repetido
    else
        # Concatenar
        cat /home/usuario/Seqs/Mepram/repositorio/01_reads/$id".fastq.gz" fastq_pass/$bc/*.fastq.gz > 01_reads/$id.fastq.gz;
        # Actualizar repo
        cp 01_reads/$id.fastq.gz /home/usuario/Seqs/Mepram/repositorio/01_reads/$id".fastq.gz"
    fi
    echo $id$'\t'$'\t'$'\t'$PWD/01_reads/$id.fastq.gz$'\t'$'\t' >> samplesheet.tsv;
    done
}
