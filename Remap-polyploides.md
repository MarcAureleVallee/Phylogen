# Ajout de variants dans les séquences consensus en fonction de la ploïdie

Ce code permet d’ajouter des variants sous forme de codes IUPAC dans les séquences consensus, tout en tenant compte de la ploïdie de l’échantillon. Avant d’exécuter ce code, il est nécessaire de disposer d’un fichier `.txt` contenant la ploïdie de chaque échantillon, sous le format : `Nom-échantillon ploïdie`

Par exemple:
```bash
Crataegus-aemula_ELB367 4
Crataegus-aestivalis_Cr454A 2
Crataegus-aestivalis_Cr477RU 2
Crataegus-ashei_ELB341 4
Crataegus-aurescens_ELB359 4
Crataegus-barroussanaVarJamensis_ELB339 4
Crataegus-berberifoliaVarBerberifolia_ELB345 2
Crataegus-brachyacantha_Cr138B 2
Crataegus-brachyacantha_Cr186TAD 2
Crataegus-brachyacantha_Cr449A 2
Crataegus-brainerdii_ELB336 4
Crataegus-calpodendron_Cr157MBG 4
Crataegus-calpodendron_Cr183TAD 2
etc...
```
Ce fichier doit être enregistré dans le répertoire de travail et être nommé `Samples-ploidy.txt`. 
```bash 
## Ajuster les variables ci-dessous de façon appropriée
SRC=/opt
WD=/scratch/$USER/HybSeqTest/remap
READS_PATH=/data/$USER/HybSeqTrim
SAMPLE_LIST_PATH=$WD/Samples-ploidy.txt
EXON=$WD/remap/70concatenated_exon_alignement
REF_FILE_ENDING=.FNA
EMAIL=votre.courriel@umontreal.ca
TIME="24-00:00:00"
CPUS=1
MEM_PER_CPU=4G

# Dossier pour remappage données Illumina sur séquences de ParalogWizard
cd $WD/remap

## Créer le batch file pour SLURM
echo '#!/bin/bash' > bwa-gatk-polyploid.sbatch
echo "#SBATCH --job-name=bwa-gatk-polyploid
#SBATCH --output=bwa-gatk-%a.out
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem-per-cpu=$MEM_PER_CPU
#SBATCH --time=$TIME

BASE_DIR=\$(pwd)

## Déterminer quel échantillon analyser dans cette tâche et définir les variables
LINE=$(head -n $SLURM_ARRAY_TASK_ID "$BASE_DIR/Samples-ploidy.txt" | tail -1)
read SAMPLE PLOIDY <<< "$LINE"


## Créer un dossier pour cet échantillon et y naviguer
mkdir -p \$BASE_DIR/\$SAMPLE
cd \$BASE_DIR/\$SAMPLE

## Créer une liste avec le nom d'échantillon (nécessaire pour seqtk)
echo \"\$SAMPLE\" > tmp.list

## Boucle qui va chercher séquence d'échantillon de chaque locus 
## et le mets dans un seul fichier .fasta, avec le nom de séquence 
## qui inclut le nom du gène
for ALIGNMENT in $EXON/*$REF_FILE_ENDING
  do
    LOCUS=\$(basename \$ALIGNMENT .fasta)
		seqtk subseq \$ALIGNMENT tmp.list >> \$SAMPLE.fasta
		sed -i \"s/\$SAMPLE.*/\$LOCUS/g\" \$SAMPLE.fasta
		sed -i \"s/\$SAMPLE multi_.*/\$LOCUS/g\" \$SAMPLE.fasta
    sed -i \"s/\$SAMPLE single_.*/\$LOCUS/g\" \$SAMPLE.fasta
  done

## Supprimer la liste avec le nom d'échantillon
rm tmp.list

## Mapper les données Illumina de l'échantillon sur les séquences références de cet échantillon
bwa index \$SAMPLE.fasta
bwa mem -t $CPUS \$SAMPLE.fasta \\
  $READS_PATH/\${SAMPLE}_trim_R1.fastq.gz \\
  $READS_PATH/\${SAMPLE}_trim_R2.fastq.gz | 
  samtools view -bS - | 
  samtools sort -o \$SAMPLE.bam -
  
## Créer les index et références nécessaires pour gatk
java -jar $SRC/picard.jar AddOrReplaceReadGroups \\
    I=\$SAMPLE.bam \\
    O=\$SAMPLE.rg.bam \\
    RGID=group1 \\
    RGLB=lib1 \\
    RGPL=illumina \\
    RGPU=unit1 \\
    RGSM=\$SAMPLE
		
samtools index \$SAMPLE.rg.bam
samtools faidx \$SAMPLE.fasta
java -jar $SRC/picard.jar CreateSequenceDictionary \\
  -R \$SAMPLE.fasta \\
  -O \$SAMPLE.dict

## Identifier les SNPs avec gatk
$SRC/gatk-4.6.0.0/gatk HaplotypeCaller \
  -R $SAMPLE.fasta \
  -I $SAMPLE.rg.bam \
  -ploidy $PLOIDY \
  -O $SAMPLE.vcf


## Filtrer les variants de mauvaise qualité
GATK_FILTER=\"DP < 3 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"
$SRC/gatk-4.6.0.0/gatk VariantFiltration \\
    -R \$SAMPLE.fasta \\
    -V \$SAMPLE.vcf \\
    -O \$SAMPLE.filt.vcf \\
    --filter-expression \"\$GATK_FILTER\" \\
    --filter-name \"StandardFilters\"

## Créer une séquence consensus avec des codes IUPAC pour les positions hétérozygotes
bgzip -c \$SAMPLE.filt.vcf > \$SAMPLE.filt.vcf.gz
tabix -p vcf \$SAMPLE.filt.vcf.gz
bcftools consensus \\
  --iupac-codes \\
  --samples \$SAMPLE \\
  --fasta-ref \$SAMPLE.fasta \\
  \$SAMPLE.filt.vcf.gz > \$SAMPLE.hetero.fasta" >> bwa-gatk.sbatch
```
## Soumettre les tâches
```bash
## Déterminer combien d'échantillons à analyser
NFILES=$(wc -l < Samples-ploidy.txt)

## Soumettre ces les tâches de mappage et identification des variants
sbatch --mail-user=$EMAIL --array=1-$NFILES bwa-gatk.sbatch
```
