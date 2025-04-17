# Extraction des données chloroplastiques à partir de données Hyb-Seq
Ce code permet d'extraire et d'assembler sur un génome de références des séquences chloroplastiques issues de données Hyb-Seq. Les lectures doivent avoir été filtrées au préalable.
## Assembler les lecures filtrées avec SPAdes
```bash
WD=/scratch/mvallee/TP_session/cp
READS=/scratch/mvallee/TP_session/data_trim
EMAIL=marcoaurelevallee@gmail.com
TIME="0-12:00:00"
MEM_PER_CPU=8G


mkdir $WD
cd $WD

echo '#!/bin/bash' > spades.sbatch
echo "#SBATCH --job-name=spades-assembly
#SBATCH --output=spades-assembly.out
#SBATCH --mail-type=END
#SBATCH --mail-user=$EMAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=$MEM_PER_CPU
#SBATCH --time=$TIME


# Sélectionner le fichier à traiter
READ1=\$(ls $READS/*trim_R1.fastq.gz | head -n \$SLURM_ARRAY_TASK_ID | tail -n 1) 

# Extraire le nom de base du fichier pour l'utiliser dans les noms de fichiers
NAME=\$(basename \$READ1 _trim_R1.fastq.gz)

# Lancer l'assemblage SPAdes avec les fichiers de lecture R1 et R2
spades.py --threads 8 --memory 60 -1 $READS/\${NAME}_trim_R1.fastq.gz -2 $READS/\${NAME}_trim_R2.fastq.gz -o $WD/\${NAME}" >> spades.sbatch
```

## Déterminer le nombre d'échantillons à analyser et envoyer la tache à slurm
```bash
NFILES=$(ls -1 $READS/*trim_R1.fastq.gz | wc -l)

# Envoyer la tache à slurm
conda activate hybpiper
sbatch --mail-user=$EMAIL --array=1-$NFILES spades.sbatch
```
# Aligner les scaffolds sur un génome de référence avec BLASTn

Les génomes de références peuvent se trouver sur le site du  [NCBI](https://www.ncbi.nlm.nih.gov/). Vous pouvez effectuer une recherche sur le site pour trouver un génome le plus près possible de votre groupe d'intéret. Ici, nous allons utiliser le génome chloroplastique de _Crataegus hupehensis_. 


## Indexer le génome de référence
Le génome de référence choisi doit avoir été au préalablement téléchargé à partir de NBCI puis envoyé dans votre environnement de travail au format .fasta. Afin d'améliorer l'assemblage, une copie de la région inversé a été supprimé dans Genious.
```bash
# Créer le dossier pour le génome de référence
mkdir /scratch/mvallee/TP_session/cp/ref_genome
cd /scratch/mvallee/TP_session/cp/ref_genome

# Indexer le génome de référence
makeblastdb -in c_hupehensis_CP-REF_whitout_IRa.fasta -dbtype nucl -out C_hupensis_cp_ref
```
## Écrire le code d'assemblage avec BLASTn
```bash
INPUT_DIR=/scratch/mvallee/TP_session/cp
OUT_DIR=/scratch/mvallee/TP_session/cp/blast
GENOME_REF=/scratch/mvallee/TP_session/cp/ref_genome/
EMAIL=marcoaurelevallee@gmail.com
TIME="0-12:00:00"
CPU=8
MEM_PER_CPU=2G

mkdir -p $OUT_DIR
cd $OUT_DIR

# Écriture du fichier SLURM
cat << EOF > blast.sbatch
#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --output=blast_%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=$EMAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$CPU
#SBATCH --mem-per-cpu=$MEM_PER_CPU
#SBATCH --time=$TIME

# Liste des dossiers
SAMPLE_DIR=\$(ls -1d $INPUT_DIR/* | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1)
SAMPLE_NAME=\$(basename \$SAMPLE_DIR)
QUERY=\$SAMPLE_DIR/scaffolds.fasta
OUTPUT=$OUT_DIR/\${SAMPLE_NAME}_blast.txt

# BLAST si le fichier existe
if [[ -f \$QUERY ]]; then
    echo "Traitement de \$SAMPLE_NAME ..."
    blastn -query \$QUERY -db $GENOME_REF -out \$OUTPUT \\
           -evalue 1e-10 -outfmt 6 -num_threads $CPU
else
    echo "Fichier manquant pour \$SAMPLE_NAME : \$QUERY"
fi
EOF
```
## Soumettre la tâche à SLURM
```bash
NFILES=$(find $INPUT_DIR -maxdepth 1 -type d -name 'Cra*' | wc -l)

conda activate hybpiper
sbatch --mail-user=$EMAIL --array=1-$NFILES blast.sbatch
```
