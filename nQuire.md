
# Estimation de la ploidie basée sur l'analyse des fréquences alléliques

## Indexer la séquence de référence
Ici, la référence utilisée est celle correspondant aux sondes Malineae481. Il faut changer la référence dépendemment des sondes utilisées lors du laboratoire (angiosperm353, Compositae-1061, etc.)
```bash
TARGETS=/home/maval/projects/def-bourret/maval/HybSeqTest/refs
cd $TARGETS

module load bwa/0.7.18
bwa index Malinae_Kew_probes_concat_exons_introns_HP.fasta
```

## Aligner les lectures sur les séquences de références 
```bash
WD="/home/maval/scratch/ploidie"
EMAIL="marcoaurelevallee@gmail.com"
TIME="0-12:00:00"
SLURM_SCRIPT="bwa-samtools.sbatch"
readspath="/home/maval/projects/def-bourret/maval/HybSeqTest/reads/trim/TEST"  
R1_files_extension=._1.trimmed.fastq.gz
R2_files_extension=._2.trimmed.fastq.gz

# Aller dans le répertoire de travail
mkdir $WD
cd $WD

# Créer le fichier de soumission SLURM
cat <<EOT > $SLURM_SCRIPT
#!/bin/bash
#SBATCH --job-name=bwa-samtools
#SBATCH --output=bwa-samtools-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=$EMAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=$TIME

# Définir les chemins
readspath="/home/maval/projects/def-bourret/maval/HybSeqTest/reads/trim/TEST"  
ploidie="$WD"      
reference="/home/maval/projects/def-bourret/maval/HybSeqTest/refs/Malinae_Kew_probes_concat_exons_introns_HP.fasta" 

# Charger les modules nécessaires
module load bwa/0.7.18
module load samtools

# Boucle sur chaque paire de fichiers contenant les lectures filtrées
for read1 in "\$readspath"/*$R1_files_extension; do
    # Identifier les fichiers R1 et R2
    read2="\${read1/$R1_files_extension/$R2_files_extension}"

    # Vérifier si le fichier R2 correspondant existe
    if [[ -f "\$read2" ]]; then
        # Extraire le nom de base du fichier sans l'extension
        base_name=\$(basename "\$read1" $R1_files_extension)

        # Créer le fichier .sam
        sam_file="\$ploidie/\${base_name}.sam"
        bwa mem "\$reference" "\$read1" "\$read2" > "\$sam_file"

        # Convertir le fichier .sam en .sorted.bam et le placer dans le répertoire ploidie
        bam_file="\$ploidie/\${base_name}.sorted.bam"
        samtools view -Sb "\$sam_file" | samtools sort -o "\$bam_file"

        # Indexer le fichier .sorted.bam
        samtools index "\$bam_file"

       # Supprimer le fichier .sam pour économiser de l'espace (optionnel)
        rm "$sam_file"

        echo "Fichier .sorted.bam créé pour \$base_name"
    else
        echo "Aucun fichier R2 trouvé pour \$read1"
    fi
done
EOT
```
## Déterminer le nombre d'échantillons à analyser et soumettre la tache à SLURM

```bash
NFILES=$(ls -1 $readspath/*$R1_files_extension | wc -l)

sbatch --mail-user=$EMAIL --array=1-$NFILES bwa-samtools.sbatch
```
## Télécharger et activer nQuire
Référez vous [ici](https://github.com/clwgg/nQuire/blob/master/README.org) pour plus d'infos.

```bash
git clone --recursive https://github.com/clwgg/nQuire
cd nQuire
make submodules
make
```
## Exécuter nQuire sur les fichiers .BAM
Assurez vous que le nquire_path correspond à l'emplacement où nQuire a été téléchargé.

```bash
cat <<"EOT" > nQuire.sbatch
#!/bin/bash
#SBATCH --job-name=nQuire
#SBATCH --output=bwa-nQuire-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=marcoaurelevallee@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-12:00:00

data_folder="/home/maval/scratch/ploidie"  
samples_file="/home/maval/scratch/ploidie/Samples_list.txt"  
nQuire_path="/home/maval/projects/def-bourret/maval/nQuire"


# Lecture du fichier de liste des échantillons
while read -r sample; do
    echo "Processing sample: $sample"
    # Utiliser le chemin complet pour nQuire
    "${nQuire_path}/nQuire" create -b "${data_folder}/${sample}.sorted.bam" -o "$sample"
    
    "${nQuire_path}/nQuire" denoise -o "${sample}_denoised" "${sample}.bin"
    
    # Déplacement des fichiers .bin vers le dossier ploidie
    mv "${sample}.bin" "${data_folder}"
    mv "${sample}_denoised.bin" "${data_folder}"
    
    echo "Finished processing $sample"
done < "$samples_file"

# Modèle de ploïdie (à exécuter après toutes les étapes ci-dessus)
denoised_bin_list=$(find "${data_folder}" -name "*_denoised.bin" -printf "%p ")
"${nQuire_path}/nQuire" lrdmodel -t 12 $denoised_bin_list > "${data_folder}/lrdmodel.tsv"
EOT
```
## Envoyer la tache à SLURM
```bash
sbatch nQuire.sbatch
```
Le résulat de ce code sera un fichier .tsv nommé `lrdmodel.tsv`. Copier le contenu de ce fichier dans excel et référez vous à [Weiß et al. (2018)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2128-z) pour interpréter les résultats.
