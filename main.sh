### 0. Inicializar carpetas y variables.

# Asegurarse de que la carpeta de trabajo tiene la carpeta fastq_pass y los archivos final_summary_*.txt, report_*.json datos_analisis.csv y datos_seq.csv

# Renombrar reads, concatenar con reads antiguos y guardar los concatenados en el repositorio:
rename_fastqs2() {
    mkdir 01_reads
    echo ID$'\t'R1$'\t'R2$'\t'LongFastQ$'\t'Fast5$'\t'GenomeSize$'\t'Fasta > samplesheet.tsv 
    tail -n +2 lista_seq.tsv | while read -r cult cep id bc c rep; do
    echo $id $bc $rep
    # Si no es repetido
    if [ "$rep" == "y" ]; then
	echo "Está repetido";
        # Concatenar
        cat ../repositorio/01_reads/$id".fastq.gz" fastq_pass/barcode$bc/*.fastq.gz > 01_reads/$id.fastq.gz;
        # Actualizar repo
        cp 01_reads/$id.fastq.gz ../repositorio/01_reads/$id".fastq.gz"
    else
        echo "No repetido, no concateno"
	# Renombrar
        cat fastq_pass/barcode$bc/*.fastq.gz > 01_reads/"$id".fastq.gz;
        # Actualizar repo
        cp 01_reads/"$id".fastq.gz ../repositorio/01_reads/$id".fastq.gz";
    fi
    echo $id$'\t'$'\t'$'\t'$PWD/01_reads/$id.fastq.gz$'\t'$'\t' >> samplesheet.tsv;
    done
    cut -f3 lista_seq.tsv | tail -n+2 > ../repositorio/samples
}

rename_fastqs2

# 0.1 Crear carpetas
mkdir 01_reads 02_filter 03_assemblies 04_taxonomies 05_plasmids 08_Anotacion 09_phages 10_ices 11_integrons

# 0.2 Meter en la lista de samples sólo los que hayan salido bien o medio bien (150Mb)
find  01_reads -size +150M | cut -d '/' -f 2 | cut -d '.' -f 1 |  sort -f | uniq > samples

# 0.3 Definir el path de las DB a usar
if [ "$ordenador" == "r" ]; then
    kraken_db="/media/usuario/datos/Databases/kraken2/"
    gambit_db="/media/usuario/datos/Databases/gambit/"
    bakta_db="/media/usuario/datos/Databases/bakta_light_5_1/db-light/"
    megares="/media/usuario/datos/Databases/megares_v3.00/megares_drugs_database_v3.00.fasta"
    isfinder_db="/media/usuario/datos/Databases/ISfinder/ISfinder-nucl.fasta"
	dorado_cmd="/home/usuario/Programs/dorado-1.1.1-linux-x64/bin/dorado"
	conda activate nano2
else
    kraken_db="/home/usuario/Databases/Kraken"
    gambit_db="/home/usuario/Databases/gambit/"
    bakta_db="/home/usuario/Databases/bakta/db/"
    megares="/home/usuario/Databases/megares/sequences"
    isfinder_db="/home/usuario/Databases/ISfinder/ISfinder-nucl.fasta"
	dorado_cmd="dorado"
fi
phastest_path="/home/usuario/Programs/phastest-docker/"


### 1. QC de lecturas.
# Nuevo entorno unificado del bloque
conda activate pitis_qc

# 1.1 Estadísticas de las lecturas crudas. Paralelizado
for i in $(cat samples); do
    NanoPlot --fastq 01_reads/${i}.fastq.gz -o 01_reads/QC/${i} --downsample 20000 --loglength --threads 3 &
done

# 1.2 Filtrar las lecturas por calidad, tamaño y eliminar inicio. Paralelizado
for i in $(cat samples); do
    gunzip   -c 01_reads/${i}.fastq.gz | chopper -q 12 -l 300 --headcrop 20 --threads 30 | gzip > 02_filter/${i}.fastq.gz &
done
wait    # Continúa con el siguiente comando solo cuando todos los procesos hayan terminado

# 1.3 Estadísticas de las lecturas filtradas. Paralelizado
for i in $(cat samples); do
    NanoPlot --fastq 02_filter/${i}.fastq.gz -o 02_filter/QC/${i} --downsample 20000 --loglength --threads 3 &
done

conda deactivate



### 2. Ensamblaje.
# Nuevo entorno unificado del bloque
conda activate pitis_ass

# 2.1 Assembly de novo
for i in $(cat samples); do
    flye --nano-hq 02_filter/${i}.fastq.gz --threads 30 --out-dir 03_assemblies/${i};
done

# 2.2 Desconcatenación de multímeros
for i in $(cat samples); do
    echo $i;
    python3 ../PITIsfinder/Scripts/deconcat/bin/deconcat.py --fasta_file 03_assemblies/${i}/assembly.fasta --fastq_file 02_filter/$i.fastq.gz --out_path 03_assemblies/$i/deconcat;
    cp 03_assemblies/${i}/deconcat/assembly_corr.fasta 03_assemblies/${i}.fasta
	rm -r 03_assemblies/${i}/deconcat/*.fastq.gz 03_assemblies/${i}/deconcat/blastn 03_assemblies/${i}/deconcat/multimer_mapping
done > 03_assemblies/deconcat.log

# 2.3 Pulir assemblies
for i in $(cat samples); do
    fq="02_filter/$i.fastq.gz"
    asm="03_assemblies/${i}/deconcat/assembly_corr.fasta"
    outdir=03_assemblies/$i/dorado_polish

    # read the first header line (robust)
    header=$(gzip -cd -- "$fq" | sed -n '1p' || true)

	# If FASTQ header indicates r9.4.1 flowcell
    if printf '%s\n' "$header" | grep -q 'dna_r9.4.1'; then
        medaka_consensus -b 4 -i 02_filter/${i}.fastq.gz -d $asm -o 03_assemblies/${i}/medaka  -m r941_e81_sup_g514 -t 30
        mv 03_assemblies/${i}/medaka/consensus.fasta 03_assemblies/${i}.fasta
        rm -r 03_assemblies/${i}/medaka/

    # If FASTQ header contains RG:Z:, run dorado with --add-fastq-rg
    elif printf '%s\n' "$header" | grep -q 'RG:Z:'; then
        echo "$i contains RG:Z: tag"
        $dorado_cmd aligner --add-fastq-rg --output-dir "$outdir" "$asm" "$fq"
        model=$(zcat $fq | head -n1 | grep -o 'dna_[^[:space:]]*' | rev | cut -f3- -d'_' | rev)
        rg=$(zcat $fq | head -n1 | cut -f3 | cut -f3 -d':')
        samtools addreplacerg -w -r "@RG\tID:${rg}\tDS:basecall_model=${model}" -o $outdir/$i.RG.bam -O bam $outdir/$i.bam
        rm $outdir/$i.bam $outdir/$i.bam.bai
        bam=$outdir/$i.RG.bam
        samtools index $bam
        $dorado_cmd polish --bacteria --batchsize 6 $bam $asm -o $outdir
        cp $outdir/consensus.fasta 03_assemblies/$i.fasta

    else
        echo "$i does NOT contain RG:Z: tag"
        $dorado_cmd aligner --output-dir "$outdir" "$asm" "$fq"
        samtools addreplacerg -r "@RG\tID:A\tDS:basecall_model=`gzip -cd "$fq" | head -n 1 | awk -F 'basecall_model_version_id=' '{print $2}'`" $outdir/$i.bam -o $outdir/$i.RG.bam -O bam
        rm $outdir/$i.bam $outdir/$i.bam.bai
        bam=$outdir/$i.RG.bam
        samtools index $bam
        $dorado_cmd polish --bacteria --batchsize 6 $bam $asm -o $outdir
        cp $outdir/consensus.fasta 03_assemblies/$i.fasta
    fi

done
rm -r .temp_dorado_model*
rm 03_assemblies/*/dorado_polish/*RG.bam*

# 2.4 Recircularizar
for i in $(cat samples); do
    circlator fixstart 03_assemblies/${i}.fasta 03_assemblies/${i}.fix && mv 03_assemblies/${i}.fix.fasta 03_assemblies/${i}.fasta  && rm 03_assemblies/*fix*
done

# 2.5 QC de assemblies
quast -o 03_assemblies/quast -t 30 --glimmer --rna-finding --conserved-genes-finding  --plots-format png 03_assemblies/*fasta
for i in $(cat samples); do
    Bandage image 03_assemblies/${i}/assembly_graph.gfa 03_assemblies/${i}.png --names --lengths --depth --fontsize 20 --toutline 3.0 --centre --query $megares --qcfilter 90 --ifilter 90 --evfilter 1e-15
done

conda deactivate


### 3. Taxonomía.
# Nuevo entorno unificado del bloque
conda activate pitis_tax

# 3.1 Taxonomía de reads con Kraken2
mkdir -p 04_taxonomies/kraken2
cp $kraken_db/*.k2d /dev/shm
for i in $(cat samples); do
    kraken2 --memory-mapping --db /dev/shm --minimum-base-quality 10 --minimum-hit-groups 100 --output 04_taxonomies/kraken2/${i}.out --use-names  --report 04_taxonomies/kraken2/${i}.report --gzip-compressed  02_filter/${i}.fastq.gz  --threads 30
done
rm /dev/shm/*.k2d

# 3.2 Asignar taxonomía con GAMBIT
gambit -d $gambit_db query -o 04_taxonomies/gambit.csv 03_assemblies/*.fasta

# 3.3 MLSTs
mlst 03_assemblies/*fasta -q > mlst.csv

# 3.4 Tipado de Klebsiellas
kleborate -a 03_assemblies/*fasta  -o 04_taxonomies/kleborate -m enterobacterales__species,klebsiella_pneumo_complex__amr,klebsiella_pneumo_complex__kaptive,klebsiella_pneumo_complex__mlst,escherichia__mlst_achtman,klebsiella_pneumo_complex__resistance_score,klebsiella_pneumo_complex__resistance_gene_count,klebsiella__ybst,klebsiella__cbst,klebsiella__abst,klebsiella__smst,klebsiella__rmst,klebsiella__rmpa2,klebsiella_pneumo_complex__virulence_score
cp 04_taxonomies/kleborate/enterobacterales__species_output.txt kleborate.tsv

# 3.5 Tipado de colis
ectyper -i 03_assemblies/*fasta -o 04_taxonomies/ectyper

conda deactivate



### 4. Anotación.
# Nuevo entorno unificado del bloque
conda activate pitis_ann

# 4.1 Anotación general con Bakta
for i in $(cat samples); do
    bakta --db $bakta_db --verbose --keep-contig-headers --output 08_Anotacion/${i} --threads 30  03_assemblies/${i}.fasta --force
done

# 4.2 Anotación de resistencias
for i in $(cat samples); do
    mkdir -p 08_Anotacion/${i}/abricate && abricate --minid 75 --mincov 75 03_assemblies/${i}.fasta > 08_Anotacion/${i}/abricate/${i}.tab
done
abricate --summary 08_Anotacion/*/abricate/*tab > 08_Anotacion/AbR.tab

conda deactivate



### 5. EGMs (aka PITIs).
# Nuevo entorno unificado del bloque
conda activate pitis_egm

# 5.1 Plásmidos
# 5.1.1 mob_recon (detección de plásmidos)
# docker pull quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1
for i in $(cat samples); do
    # docker run --rm -v $(pwd):/mnt/ "kbessonov/mob_suite:3.0.3" mob_recon -i /mnt/03_assemblies/${i}.fasta -o /mnt/08_Anotacion/${i}/mob_recon -c --force -t -n 30;
    docker run --rm -u $(id -u):$(id -g) -v $(pwd):/mnt/ quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1 mob_recon -i /mnt/03_assemblies/${i}.fasta -o /mnt/08_Anotacion/${i}/mob_recon -c --force -n 30;
done  # Tipado de Plásmidos

# 5.1.2 Copla (Asignar PTU)
echo "Copla results" > copla.txt
for j in $(cat samples); do
    for i in $(find 08_Anotacion/${j}/mob_recon/*fasta  -size -600k -size +1k); do # Aplicar copla a los contigs pequeños (<600kb), sospechosos de ser plásmidos
        new_name=$(echo $i | rev | cut -f1 -d'/' | cut -f2- -d'.' | cut -f1 -d'_' | rev)    # Sacar el grupo del plásmido
        cp ${i} 05_plasmids/${new_name}_${j}.fasta
		echo "Sample: ${j}" && echo "Contig: ${i:(-11):5}" &&  docker run --rm -v $(pwd):/tmp rpalcab/copla:1.0 copla /tmp/${i} /data/app/databases/Copla_RS84/RS84f_sHSBM.pickle /data/app/databases/Copla_RS84/CoplaDB.fofn /tmp/08_Anotacion/${j}/copla
    done  >> copla.txt
done

# 5.2 Integrones
for i in $(cat samples); do
    integron_finder 03_assemblies/${i}.fasta --cpu 30 --outdir 11_integrons/${i} --func-annot --gbk;
done
python3 ../PITIsfinder/Scripts/integrones/integron_parser.py .
cp 11_integrons/integron_summary.csv .

# 5.3 Fagos
for i in $(cat samples); do
    mkdir -p 09_phages/phastest_deep/"$i"/;
    cp 03_assemblies/"$i".fasta $phastest_path/phastest_inputs/;
    docker compose -f $phastest_path/docker-compose.yml down --remove-orphans;
    docker compose -f $phastest_path/docker-compose.yml run phastest -i fasta -m deep -s "$i".fasta --phage-only --yes;
    cp -r $phastest_path/phastest-app-docker/JOBS/"$i" 09_phages/phastest_deep/;
    # rm -r $phastest_path/phastest-app-docker/JOBS/"$i";
    rm -r $phastest_path/phastest_inputs/"$i".fasta;
done

python3 ../PITIsfinder/Scripts/phages/phage_parser.py .
cp 09_phages/phage_summary.csv .

# 5.4 IS
for i in $(cat samples) ; do
    sed -i 's/>contig/>Chr/g' 08_Anotacion/${i}/mob_recon/chromosome.fasta &&
    makeblastdb -in 08_Anotacion/${i}/mob_recon/chromosome.fasta -dbtype nucl &&
    blastn -db  08_Anotacion/${i}/mob_recon/chromosome.fasta -query $isfinder_db -outfmt "6 qseqid sseqid sstart send pident mismatch evalue" | sort -k 3 > 08_Anotacion/${i}/IS_chr.tsv &&
    python3 ../PITIsfinder/Scripts/IS_parser.py -i 08_Anotacion/${i};
done

conda deactivate



### 6. Informes

# 6.1 QC de reads
# 6.1.1 QC de reads crudos
for i in $(cat samples); do
    grep 'Median read length\|Median read quality\|Number of reads\|Total bases\|1:' 01_reads/QC/${i}/NanoStats.txt | cut  -d ':' -f 2 | cut -d '(' -f 1 | column -t |datamash transpose
done > 01_reads/QC_reads.csv
paste samples 01_reads/QC_reads.csv > 01_reads/QC_reads_pre.csv
sed -i $'1 i\\\nSample\tMedian\ length\tMedian\ quality\tTotal\ reads\tTotal\ bases\tMaxQ\tLongest\ read' 01_reads/QC_reads_pre.csv
cp 01_reads/QC_reads_pre.csv QC_reads_pre.csv

# 6.1.2 QC de reads postfiltrado
for i in $(cat samples); do
    grep 'Median read length\|Median read quality\|Number of reads\|Total bases\|1:' 02_filter/QC/${i}/NanoStats.txt | cut  -d ':' -f 2 | cut -d '(' -f 1 | column -t |datamash transpose
done > 02_filter/QC_reads.csv
paste samples 02_filter/QC_reads.csv > 02_filter/QC_reads_post.csv
sed -i $'1 i\\\nSample\tMedian\ length\tMedian\ quality\tTotal\ reads\tTotal\ bases\tMaxQ\tLongest\ read' 02_filter/QC_reads_post.csv
cp 02_filter/QC_reads_post.csv QC_reads_post.csv

# 6.1.3 QC de reads unificado
paste QC_reads_pre.csv QC_reads_post.csv > QC_reads.csv
rm QC_reads_pre.csv QC_reads_post.csv


# 6.2 Informe de Taxonomía
# 6.2.1 Kraken2
for i in $(cat samples) ; do    # coge los géneros con más del 20% de reads o especies con más del 5%
    awk '($1 >= 20 && $4 == "G")'  04_taxonomies/kraken2/${i}.report | cut -f 1,6 | tr -s '  '| awk '{print i,"\t",$1,"\t",$2,$3,$4}' "i=${i}" >> 04_taxonomies/kraken2/genus.csv
    awk '($1 >= 4  && $4 == "S")'  04_taxonomies/kraken2/${i}.report | cut -f 1,6 | tr -s '  '| awk '{print i,"\t",$1,"\t",$2,$3,$4}' "i=${i}" >> 04_taxonomies/kraken2/species.csv
    awk '($1 >= 4  && $4 == "S2")' 04_taxonomies/kraken2/${i}.report | cut -f 1,6 | tr -s '  '| awk '{print i,"\t",$1,"\t",$2,$3,$4}' "i=${i}" >> 04_taxonomies/kraken2/strains.csv
done
# mergear tablas
awk -F"\t" 'NR==FNR{a[$1]=$0;next} ($1 in a){b=$1;$1="";print a[b] $0}' OFS="\t" 04_taxonomies/kraken2/genus.csv 04_taxonomies/kraken2/species.csv > 04_taxonomies/kraken_report.csv
# introducir header
sed -i $'1 i\\\nSample\tReads\(\%\)\tGenus\tReads\(\%\)\tSpecies' 04_taxonomies/kraken_report.csv
# eliminar espacios en algunas columnas y copiar en inicio
sed  's/\t[[:blank:]]*/\t/g' 04_taxonomies/kraken_report.csv > kraken.csv

# 6.3 QC de Assemblies
cut 03_assemblies/quast/transposed_report.tsv -f 1,14-17,26 > QC_assembly.csv
sed -i 's/Assembly/Samples/' QC_assembly.csv

# 6.4 Informe de Abricate
cp 08_Anotacion/AbR.tab AbR_report.csv

# 6.5 Resto del informe
python3 ../PITIsfinder/Scripts/parser.py -i .

# 6.6 Informe IS
for i in $(cat samples); do cat 08_Anotacion/${i}/IS_chr_out.tsv | cut -f 2 | sort | uniq -c | sort -r | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]/\t&/g'  >08_Anotacion/${i}/N_IS_${i}.tsv; done &&
for i in $(cat samples); do cat 08_Anotacion/${i}/N_IS* | head -n 1  >> 08_Anotacion/max_N_IS.tsv ; done && paste samples 08_Anotacion/max_N_IS.tsv > max_IS.tsv &&
for i in $(cat samples); do tail 08_Anotacion/${i}/IS_chr_out.tsv -n +2 | wc -l >> 08_Anotacion/total_N_IS.tsv ; done && paste max_IS.tsv 08_Anotacion/total_N_IS.tsv > IS.tsv &&
sed  -i '1i sample      max     IS_name Total_IS' IS.tsv &&
rm max_IS.tsv 08_Anotacion/total_N_IS.tsv 08_Anotacion/max_N_IS.tsv

# Python para la generación de datos_seq_nuevos.csv y datos_analisis_nuevos.csv
python3 ../PITIsfinder/Scripts/Datos_seq_unified2.py --input_path .

###############
# Repositorio #
###############
# Los fastqs de todas las muestras se han añadido al principio del lablog
# Copiar assemblies de las muestras que han salido bien
for i in $(cat samples); do
    # cp -r 03_assemblies/"$i".fasta ../repositorio/03_assemblies/"$i".fasta;   # solo assembly
    cp -r 03_assemblies/"$i"* ../repositorio/03_assemblies/;                    # assembly y subcarpeta
done

# Copiar anotaciones de las muestras que han salido bien (anotación general y AbR)
for i in $(cat samples); do
    mkdir ../repositorio/08_Anotacion/$i
    ls 08_Anotacion/$i/*.* | grep -v -E 'IS|resumen' | xargs cp -t ../repositorio/08_Anotacion/$i
    cp -r 08_Anotacion/$i/abricate ../repositorio/08_Anotacion/$i
done

# Copiar PITIS de las muestras que han salido bien
# Plásmidos (queremos copiar la info de MOBsuite y de copla?)
for i in $(cat samples); do
    rm -r ../repositorio/05_plasmids/*"$i"*
done
cp -r 05_plasmids/* ../repositorio/05_plasmids

# Integrones
for i in $(cat samples); do
    rm -r ../repositorio/11_integrons/*"$i"*
    cp -r 11_integrons/* ../repositorio/11_integrons/
done

# Fagos
for i in $(cat samples); do
    rm -r ../repositorio/09_phages/*"$i"*
done
cp -r 09_phages/* ../repositorio/09_phages/
