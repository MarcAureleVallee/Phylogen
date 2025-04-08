# Assembler les lecures filtrées avec SPAdes
Ce code permet d'extraire et d'assembler les séquences chloroplastiques à partir de données Hyb-Seq. Les lectures doivent avoir été filtrées au préalable.
```bash
WD=/scratch/mvallee/TP_session/cp
READS=/scratch/mvallee/TP_session/data_trim
EMAIL=marcoaurelevallee@gmail.com
TIME="0-12:00:00"
MEM_PER_CPU=2G


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
spades.py -1 $READS/\${NAME}_trim_R1.fastq.gz -2 $READS/\${NAME}_trim_R2.fastq.gz --isolate -o $WD" >> spades.sbatch
```

## Déterminer le nombre d'échantillons à analyser et envoyer la tache à slurm
```bash
NFILES=$(ls -1 $READS/*trim_R1.fastq.gz | wc -l)

# Envoyer la tache à slurm
conda activate hybpiper
sbatch --mail-user=$EMAIL --array=1-$NFILES spades.sbatch
```
