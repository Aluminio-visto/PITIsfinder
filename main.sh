cp -r /var/lib/minknow/data/$CARPETA_RUN/no_sample_id/*/fastq_pass/ .

## Copiar final_summary_*.txt y report_*.json a la carpeta de trabajo
cp /var/lib/minknow/data/$CARPETA_RUN/no_sample_id/*/final_summary_*.txt /var/lib/minknow/data/CARPETA_RUN/no_sample_id/*/report_*.json .

## Copiar datos_analisis_nuevos.csv y datos_seq_nuevos.csv a la carpeta actual, quitándoles el sufijo _nuevos
last_run=$(ls -d ../202*/ | tail -n2 | head -n1)
cp $last_run/datos_analisis_nuevos.csv datos_analisis.csv
cp $last_run/datos_seq_nuevos.csv datos_seq.csv

# Renombrar reads, concatenar con reads antiguos y guardar los concatenados en el repositorio:
rename_fastqs() {
    mkdir 01_reads
    echo ID$'\t'R1$'\t'R2$'\t'LongFastQ$'\t'Fast5$'\t'GenomeSize$'\t'Fasta > samplesheet.tsv
    tail -n +2 lista_seq.tsv | while read -r cult cep id bc c rep; do
    echo $id $bc $rep
    # Si no es repetido
    if [ -z "$rep" ]; then
        # Renombrar
        cat fastq_pass/barcode$bc/*.fastq.gz > 01_reads/"$id".fastq.gz;
        # Actualizar repo
        cp 01_reads/"$id".fastq.gz ../repositorio/01_reads/$id".fastq.gz";
    # Si es repetido
    else
        # Concatenar
        cat ../repositorio/01_reads/$id".fastq.gz" fastq_pass/barcode$bc/*.fastq.gz > 01_reads/$id.fastq.gz;
        # Actualizar repo
        cp 01_reads/$id.fastq.gz ../repositorio/01_reads/$id".fastq.gz"
    fi
    echo $id$'\t'$'\t'$'\t'$PWD/01_reads/$id.fastq.gz$'\t'$'\t' >> samplesheet.tsv;
    done
    cut -f3 lista_seq.tsv | tail -n+2 > ../repositorio/samples
}

rename_fastqs

# Crear carpetas
mkdir 01_reads 02_filter 03_assemblies 04_taxonomies 05_plasmids 08_Anotacion 09_phages 10_ices 11_integrons

# Meter en la lista de samples sólo los que hayan salido bien o medio bien (10Mb)
find  01_reads -size +150M | cut -d '/' -f 2 | cut -d '.' -f 1 |  sort -f | uniq > samples

if [ "$ordenador" == "r" ]; then
    conda activate nano2
    kraken_db="/media/usuario/datos/Databases/kraken2/"
    gambit_db="/media/usuario/datos/Databases/gambit/"
    bakta_db="/media/usuario/datos/Databases/bakta/db/"
    megares="/media/usuario/datos/Databases/megares/sequences"
else
    conda activate minion
    kraken_db="/home/usuario/Databases/Kraken"
    gambit_db="/home/usuario/Databases/gambit/"
    bakta_db="/home/usuario/Databases/bakta/db/"
    megares="/home/usuario/Databases/megares/sequences"
fi

# QC
for i in $(cat samples); do
    NanoPlot --fastq 01_reads/${i}.fastq.gz -o 01_reads/QC/${i} --downsample 20000 --loglength --threads 3 &
done

# Filtrar
for i in $(cat samples); do
    gunzip   -c 01_reads/${i}.fastq.gz | chopper -q 12 -l 300 --headcrop 20 --threads 30 | gzip > 02_filter/${i}.fastq.gz &
done
wait    # Continúa con el siguiente comando solo cuando todos los procesos hayan terminado

# QC
for i in $(cat samples); do
    NanoPlot --fastq 02_filter/${i}.fastq.gz -o 02_filter/QC/${i} --downsample 20000 --loglength --threads 3 &
done

# Taxonomía de reads
if [ "$ordenador" = "r" ]; then
    # Tristemente, en mi ordenador no se puede cargar la db en memoria, peta.
    for i in $(cat samples); do
        kraken2 --db $kraken_db --minimum-base-quality 10 --minimum-hit-groups 100 --output 04_taxonomies/kraken2/${i}.out --use-names  --report 04_taxonomies/kraken2/${i}.report --gzip-compressed  02_filter/${i}.fastq.gz  --threads 30
    done
else
    mkdir -p 04_taxonomies/kraken2
    cp $kraken_db/*.k2d /dev/shm
    for i in $(cat samples); do
        kraken2 --memory-mapping --db /dev/shm --minimum-base-quality 10 --minimum-hit-groups 100 --output 04_taxonomies/kraken2/${i}.out --use-names  --report 04_taxonomies/kraken2/${i}.report --gzip-compressed  02_filter/${i}.fastq.gz  --threads 30
    done
    rm /dev/shm/*.k2d
fi

# Assembly rápido
for i in $(cat samples); do
    flye --nano-hq 02_filter/${i}.fastq.gz --threads 30 --out-dir 03_assemblies/${i};
done

# Desconcatenación de multímeros
conda activate deconcat
for i in $(cat samples); do
    echo $i;
    python3 ../Scripts/deconcat/deconcat.py --fasta_file 03_assemblies/${i}/assembly.fasta --fastq_file 02_filter/$i.fastq.gz --out_path 03_assemblies/$i/deconcat;
    cp 03_assemblies/${i}/deconcat/assembly_corr.fasta 03_assemblies/${i}.fasta
done > 03_assemblies/deconcat.log
conda deactivate

# Pulir assemblies rápidos
conda activate  medaka
for i in $(cat samples); do
    medaka_consensus -b 4 -i 02_filter/${i}.fastq.gz -d 03_assemblies/${i}.fasta -o 03_assemblies/${i}/medaka  -m r941_e81_sup_g514 -t 30
    mv 03_assemblies/${i}/medaka/consensus.fasta 03_assemblies/${i}.fasta
    rm -r 03_assemblies/${i}/medaka/
done
conda deactivate

# Recircularizar
if [ "$ordenador" = "r" ]; then
    :
else
    conda activate ~/miniconda3/envs/nano
fi

for i in $(cat samples); do
    circlator fixstart 03_assemblies/${i}.fasta 03_assemblies/${i}.fix && mv 03_assemblies/${i}.fix.fasta 03_assemblies/${i}.fasta  && rm 03_assemblies/*fix*
done
if [ "$ordenador" = "r" ]; then
    :
else
    conda deactivate
fi

#QC de assemblies
/home/usuario/Programs/quast-5.2.0/quast.py -o 03_assemblies/quast -t 30 --glimmer --rna-finding --conserved-genes-finding  --plots-format png 03_assemblies/*fasta
for i in $(cat samples); do
    Bandage image 03_assemblies/${i}/assembly_graph.gfa 03_assemblies/${i}.png --names --lengths --depth --fontsize 20 --toutline 3.0 --centre --query $megares --qcfilter 90 --ifilter 90 --evfilter 1e-15
done

# Asignar taxonomía con GAMBIT
if [ "$ordenador" = "r" ]; then
    conda activate gambit
else
    :
fi
gambit -d $gambit_db query -o 04_taxonomies/gambit.csv 03_assemblies/*.fasta
mkdir -p 04_taxonomies/gtdb
cp 03_assemblies/*fasta 04_taxonomies/gtdb/

if [ "$ordenador" = "r" ]; then
    conda deactivate
    conda activate icecreen
else
    :
fi

# Anotación con Bakta
for i in $(cat samples); do
    bakta --db $bakta_db --verbose --output 08_Anotacion/${i} --threads 30  03_assemblies/${i}.fasta --force
done

# Sacar plásmidos con mob_recon e integrones con integron_finder
if [ "$ordenador" = "r" ]; then
    conda deactivate
else
    conda activate ~/miniconda3/envs/covidion
fi
for i in $(cat samples); do
    docker run --rm -v $(pwd):/mnt/ "kbessonov/mob_suite:3.0.3" mob_recon -i /mnt/03_assemblies/${i}.fasta -t -o /mnt/08_Anotacion/${i}/mob_recon -c --force -t -n 30;
done  # Tipado de Plásmidos

for i in $(cat samples); do
    integron_finder 03_assemblies/${i}.fasta --cpu 30 --outdir 11_integrons/${i} --func-annot --gbk;
done
if [ "$ordenador" = "r" ]; then
    :
else
    conda deactivate
fi
python3 ../Scripts/integrones/integron_parser.py .
cp 11_integrons/integron_summary.csv .

# Aplicar copla a los contigs pequeños (<600kb), sospechosos de ser plásmidos
conda activate copla
echo "Copla results" > copla.txt
for j in $(cat samples); do
    for i in $(find 08_Anotacion/${j}/mob_recon/*fasta  -size -600k -size +1k); do
        new_name=$(echo $i | rev | cut -f1 -d'/' | cut -f2- -d'.' | cut -f1 -d'_' | rev)    # Sacar el grupo del plásmido
        cp ${i} 05_plasmids/${new_name}_${j}.fasta
        echo "Sample: ${j}" && echo "Contig: ${i:(-11):5}" &&  python3 ~/Programs/copla/bin/copla.py ${i} /home/usuario/Programs/copla/databases/Copla_RS84/RS84f_sHSBM.pickle /home/usuario/Programs/copla/databases/Copla_RS84/CoplaDB.fofn 08_Anotacion/${j}/copla
    done  >> copla.txt
done
conda deactivate
# Renombrar contigs en función de si son cromosoma o plásmido, añadiendo nombre de muestra

# Sacar resistencias
conda activate abr
for i in $(cat samples); do
    mkdir -p 08_Anotacion/${i}/abricate && abricate --minid 75 --mincov 75 03_assemblies/${i}.fasta > 08_Anotacion/${i}/abricate/${i}.tab
done
abricate --summary 08_Anotacion/*/abricate/*tab > 08_Anotacion/AbR.tab
conda deactivate

# Sacar MLSTs
if [ "$ordenador" = "r" ]; then
    conda activate mlst
else
    conda activate ~/miniconda3/envs/mlst
fi
mlst 03_assemblies/*fasta -q > mlst.csv
conda deactivate

# Tipar bichos particulares
# Klebsiellas y colis
conda activate kleb
kleborate -a 03_assemblies/*fasta  -o klebsiellas -m enterobacterales__species,klebsiella_pneumo_complex__amr,klebsiella_pneumo_complex__kaptive,klebsiella_pneumo_complex__mlst,escherichia__mlst_achtman,klebsiella_pneumo_complex__resistance_score,klebsiella_pneumo_complex__resistance_gene_count,klebsiella__ybst,klebsiella__cbst,klebsiella__abst,klebsiella__smst,klebsiella__rmst,klebsiella__rmpa2,klebsiella_pneumo_complex__virulence_score
cp klebsiellas/enterobacterales__species_output.txt kleborate.tsv
ectyper -i 04_taxonomies/gtdb -o 04_taxonomies/ectyper
conda deactivate

# Fagos
for i in $(cat samples); do
    mkdir -p 09_phages/phastest_deep/"$i"/;
    sudo cp 03_assemblies/"$i".fasta /home/usuario/Programs/phastest-docker/phastest_inputs/;
    docker compose -f /home/usuario/Programs/phastest-docker/docker-compose.yml down --remove-orphans;
    docker compose -f /home/usuario/Programs/phastest-docker/docker-compose.yml run phastest -i fasta -m deep -s "$i".fasta --phage-only --yes;
    sudo cp -r /home/usuario/Programs/phastest-docker/phastest-app-docker/JOBS/"$i" 09_phages/phastest_deep/;
    sudo rm -r /home/usuario/Programs/phastest-docker/phastest-app-docker/JOBS/"$i";
    sudo rm -r /home/usuario/Programs/phastest-docker/phastest_inputs/"$i".fasta;
done

python3 ../Scripts/phages/phage_parser.py .
cp 09_phages/phage_summary.csv .

# # ICEs
# conda activate icecreen
# # Copiamos los cromosomas a la carpeta 10_ices
# for i in $(cat samples); do
#     cp 08_Anotacion/$i/mob_recon/chromosome.fasta 10_ices/$i".fasta"
# done

# cd 10_ices
# ls *.fasta > samples.txt
# bash ../../Scripts/ices/ice_characterize.sh -i samples.txt -o anotacion
# python3 ../../Scripts/ices/ice_parser.py -p ./
# conda deactivate
# cd ..
# cp 10_ices/tables/ICE_summary.csv .


##################################
# Generar tablas para el informe #
##################################

### QC de Reads ###
# QC de reads prefiltrado
for i in $(cat samples); do
    grep 'Median read length\|Median read quality\|Number of reads\|Total bases\|1:' 01_reads/QC/${i}/NanoStats.txt | cut  -d ':' -f 2 | cut -d '(' -f 1 | column -t |datamash transpose
done > 01_reads/QC_reads.csv
paste samples 01_reads/QC_reads.csv > 01_reads/QC_reads_pre.csv
sed -i $'1 i\\\nSample\tMedian\ length\tMedian\ quality\tTotal\ reads\tTotal\ bases\tMaxQ\tLongest\ read' 01_reads/QC_reads_pre.csv
cp 01_reads/QC_reads_pre.csv QC_reads_pre.csv

# QC de reads postfiltrado
for i in $(cat samples); do
    grep 'Median read length\|Median read quality\|Number of reads\|Total bases\|1:' 02_filter/QC/${i}/NanoStats.txt | cut  -d ':' -f 2 | cut -d '(' -f 1 | column -t |datamash transpose
done > 02_filter/QC_reads.csv
paste samples 02_filter/QC_reads.csv > 02_filter/QC_reads_post.csv
sed -i $'1 i\\\nSample\tMedian\ length\tMedian\ quality\tTotal\ reads\tTotal\ bases\tMaxQ\tLongest\ read' 02_filter/QC_reads_post.csv
cp 02_filter/QC_reads_post.csv QC_reads_post.csv

paste QC_reads_pre.csv QC_reads_post.csv > QC_reads.csv
rm QC_reads_pre.csv QC_reads_post.csv

### Informe de Kraken2 ###
# coge los géneros con más del 20% de reads o especies con más del 5%
for i in $(cat samples) ; do
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

### QC de Assemblies ###
sed '1i Samples' samples > samples2
paste samples2  03_assemblies/quast/transposed_report.tsv | cut -f 1,15-18,27 > QC_assembly.csv

### Informe de Abricate ###
cp 08_Anotacion/AbR.tab AbR_report.csv

# MLST
# Para sacar los alelos que dan fallo (y solo los que dan fallo)
python3 ../Scripts/parser.py -i .

# TxSS
for i in $(cat samples); do
    { echo ${i}_genes ;tail +6 08_Anotacion/${i}/TXSS/best_solution.tsv | head -n -2 | cut -f 3 ;} > 08_Anotacion/${i}/resumen_TXSS && paste 08_Anotacion/*/resumen_TXSS  > TXSS_resumen.tsv
done
# faltaría hacer que se organizasen por grupos con PANDAS de tal forma que quedasen a la misma altura los distintos sistemas de secreción.
# podría hacerlo generando un dataframe por cada txss y luego concatenandolos todos uno encima del otro.

# IS
for i in $(cat samples) ; do
sudo sed -i 's/>contig/>Chr/'g 08_Anotacion/${i}/mob_recon/chromosome.fasta &&
sudo makeblastdb -in 08_Anotacion/${i}/mob_recon/chromosome.fasta -dbtype nucl &&
blastn -db  08_Anotacion/${i}/mob_recon/chromosome.fasta -query /home/usuario/Databases/ISfinder/ISfinder-nucl.fasta  -outfmt "6 qseqid sseqid sstart send pident mismatch evalue" | sort -k 3 > 08_Anotacion/${i}/IS_chr.tsv &&
python3 /home/usuario/Seqs/Servicio/IS_parser.py -i 08_Anotacion/${i};
done
for i in $(cat samples); do cat 08_Anotacion/${i}/IS_chr_out.tsv | cut -f 2 | sort | uniq -c | sort -r | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]/\t&/g'  >08_Anotacion/${i}/N_IS_${i}.tsv; done &&
for i in $(cat samples); do cat 08_Anotacion/${i}/N_IS* | head -n 1  >> 08_Anotacion/max_N_IS.tsv ; done && paste samples 08_Anotacion/max_N_IS.tsv > max_IS.tsv &&
for i in $(cat samples); do tail 08_Anotacion/${i}/IS_chr_out.tsv -n +2 | wc -l >> 08_Anotacion/total_N_IS.tsv ; done && paste max_IS.tsv 08_Anotacion/total_N_IS.tsv > IS.tsv &&
sed  -i '1i sample      max     IS_name Total_IS' IS.tsv &&
rm max_IS.tsv 08_Anotacion/total_N_IS.tsv 08_Anotacion/max_N_IS.tsv

# Python para la generación de datos_seq_nuevos.csv y datos_analisis_nuevos.csv
python3 ../Scripts/Datos_seq_unified2.py --input_path .

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
cp 05_plasmids/ ../repositorio/05_plasmids

# # ICEs
# for i in $(cat samples); do
#     rm -r ../repositorio/10_ices/*/*"$i"*
#     rm -r ../repositorio/10_ices/*"$i"*
# done
# ls -d 10_ices/*/ | grep -v "anotacion" | xargs cp -r -t ../repositorio/10_ices

# Integrones
for i in $(cat samples); do
    rm -r ../repositorio/11_integrons/*"$i"*
    cp -r 11_integrons/* ../repositorio/11_integrons/
done

# Fagos
for i in $(cat samples); do
    rm -r ../repositorio/09_phages/*"$i"*
done
cp -r 09_phages/ ../repositorio/09_phages/
