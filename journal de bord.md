# Créer des alignements orthologues avec ParalogWizard
Cette étape est effectuée sur le cluster Metacentrum

## S'assurer que tous les fichiers dédupliqués et filtrés sont dans le dossier 10deduplicated reads
Voici la liste d'échantillons avec lesquelles je veux travailler

Crataegus-aemula_HbPr33
Crataegus-aestivalis_Cr477RU
Crataegus-ashei_HbPr07
Crataegus-aurescens_HbPr25
Crataegus-barroussanaVarJamensis_HbPr05
Crataegus-berberifoliaVarBerberifolia_HbPr11
Crataegus-brachyacantha_Cr449A
Crataegus-brainerdii_FtPr02
Crataegus-calpodendron_Cr462RU
Crataegus-castlegarensis_HbPr14
Crataegus-chlorosarca_Cr196LE
Crataegus-chrysocarpaVarChrysocarpa_FtPr40
Crataegus-chrysocarpaVarPraecox_FtPr41
Crataegus-chrysocarpaVarSubrotundifolia_FtPr42
Crataegus-coccineaVarCoccinea_FtPr09
Crataegus-coccineaVarPringlei_FtPr10
Crataegus-crus-galli_Cr481RU
Crataegus-cuprina_HbPr06
Crataegus-dahurica_Cr743LE
Crataegus-dodgei_FtPr43
Crataegus-douglassii_HbPr15
Crataegus-enderbyensis_HbPr39
Crataegus-flabellataVarGrayana_FtPr47
Crataegus-fluviatilis_FtPr46
Crataegus-frugiferens_HbPr03
Crataegus-furtiva_HbPr21
Crataegus-germanica_Cr001HBV
Crataegus-gracilior_HbPr12
Crataegus-grandifoliaVarGrandifolia_HbPr16
Crataegus-greggianaVarGreggiana_HbPr17
Crataegus-harbisonii_HbPr08
Crataegus-hupehensis_Cr536MO
Crataegus-intricataVarIntricata_HbPr19
Crataegus-irrasa_FtPr44
Crataegus-jonesiae_HbPr37
Crataegus-lacrimata_HbPr22
Crataegus-macracantha_FtPr23
Crataegus-macrosperma_FtPr48
Crataegus-marshalii_Cr487RU
Crataegus-mendosa_HbPr35
Crataegus-mexicana_Cr517MO
Crataegus-mollis_Cr473RU
Crataegus-monogyna_Cr175BRNU
Crataegus-nananixonii_HbPr20
Crataegus-nigra_Cr089BGBM
Crataegus-okanaganensis_HbPr38
Crataegus-opaca_Cr466RU
Crataegus-pennsylvanica_FtPr28
Crataegus-pentagyna_Cr195RU
Crataegus-phaenopyrum_Cr495RU
Crataegus-pinnatifida_Cr100BGBM
Crataegus-populnea_HbPr32
Crataegus-pulcherrima_Cr502RU
Crataegus-punctata_Cr496RU
Crataegus-purpurella_HbPr29
Crataegus-pycnoloba_Cr290RU
Crataegus-rivulodamensis_HbPr30
Crataegus-roseiSubspRoseiVarAmoena_HbPr13
Crataegus-rufula_HbPr01
Crataegus-saligna_Cr559NY
Crataegus-sanguinea_Cr813LE
Crataegus-schuettei_Cr499RU
Crataegus-serratissima_HbPr18
Crataegus-sororia_HbPr04
Crataegus-spathulata_Cr488RU
Crataegus-submollis_FtPr27
Crataegus-suborbiculata_FtPr34
Crataegus-succulentaVarSucculenta_FtPr24
Crataegus-sulfurea_HbPr45
Crataegus-texana_Cr492RU
Crataegus-tracyi_Cr490RU
Crataegus-tracyiVarMadrensis_HbPr26
Crataegus-triflora_Cr459RU
Crataegus-uniflora_HbPr31
Crataegus-venustaVarVenusta_HbPr36
Crataegus-viridis_Cr235TAD
Hesperomeles_obtusifolia_SRR12879568
Malus-domestica_Mlu855Rosi
Photinia-bodinieri_Pht724MO


## Assemblage
```bash
./ParalogWizard_1a_CastSubmitAssemble_subitter.sh
```
## Retrieving
 ```bash
qsub ParalogWizard_1b_CastSubmitRetrieve.sh
```

## Detection des paralogues
```bash
qsub ParalogWizard_2a_CastSubmitAnalyze.sh
```

## Créer des alignements orthologues
```bash
##Copier la référence personnalisée dans le répertoire source
cp datadir/41detected_par/customized_reference_div_4.8_10.8.fas HybSeqSource

qsub ParalogWizard_3_CastSubmitSeparate.sh
```
Une fois les paralogues séparer, télécharger les séquences sur l'ordinateur local puis sur le cluster amphidzen
```bash
rsync --progress marc-aurele@tarkil.metacentrum.cz:/storage/brno2/home/marc-aurele/datadir/70concatenated_exon_alignments/*.fasta mnt/c/Users/Utilisateur/Documents/SCRATCH

rsync --progress mnt/c/Users/Utilisateur/Documents/SCRATCH/*.fasta mvallee@aphidzen.irbv.umontreal.ca:/scratch/mvallee/TP_session/Assem_avec_para/exon
```

## Reformater le nom des gènes
```bash
cd /scratch/mvallee/TP_session/Assem_avec_para/exon

for file in Assembly_*.fasta; do
    mv "$file" "${file#Assembly_}"
done

##Changer .fasta en .FNA
for file in *.fasta; do
    mv "$file" "${file%.fasta}.FNA"
done

##Ajouter single_hit a la fin des noms de séquence pour que mes données soient compatibles avec le code d'ELB

for FILE in /scratch/mvallee/TP_session/Assem_avec_para/exon/*.FNA; do
    sed -i 's/^>\(.*\)/>\1 single_hit/' "$FILE"
done
for file in *.fasta; do
    mv "$file" "${file%.fasta}.FNA"
done
```
## Ajouter single_hit a la fin des noms de séquence pour que mes données soient compatibles avec le code d'ELC
```bash
for FILE in /scratch/mvallee/TP_session/Assem_avec_para/exon/*.FNA; do
    sed -i 's/^>\(.*\)/>\1 single_hit/' "$FILE"
done
```
# Remapper les lectures sur les exons pour obtenir des séquences avec codes IUPAC

## Créer un liste des échantillons a analyser
Copier dans ce fichier la liste des échantillons plus haut
```bash
nano sampleslist.txt
```

## Remplacer les caractères indéterminés par des N
```bash
input_dir="/scratch/mvallee/TP_session/Assem_avec_para/exon"
output_dir="/scratch/mvallee/TP_session/Assem_avec_para/exon_clean"

mkdir -p "$output_dir"

# Boucle pour traiter chaque fichier dans le répertoire `exon`
for file in "$input_dir"/*.FNA; do
    # Chemin complet pour le fichier nettoyé dans le répertoire cible
    cleaned_file="$output_dir/$(basename "$file")"
    
    # Remplacer les caractères "?" et "-" par "N" dans les séquences uniquement (ne modifie pas les en-têtes)
    sed -e '/^>/!s/[?-]/N/g' "$file" > "$cleaned_file"
    
    echo "Transformation terminée pour : $(basename "$file") -> $(basename "$cleaned_file")"
done
```
## Lancer BWA et gatk
```bash
# Ajuster les variables ci-dessous de façon appropriée
SRC=/opt
WD=/scratch/mvallee/TP_session/Assem_avec_para
READS_PATH=/scratch/mvallee/TP_session/data_trim
REF_TYPE=exon_clean
REF_FILE_ENDING=.FNA
EMAIL=marcoaurelevallee@gmail.com
TIME="0-12:00:00"

# Créer le dossier pour le remappage des données Illumina sur les séquences HybPiper
mkdir -p $WD/remap
cd $WD/remap

## Créer un lien vers la liste des échantillons à analyser
ln -s $WD/samplelist.txt .

## Créer le batch file pour SLURM
echo '#!/bin/bash' > bwa-gatk.sbatch
echo "#SBATCH --job-name=bwa-gatk
#SBATCH --output=bwa-gatk-%a.out
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=$TIME

BASE_DIR=\$(pwd)

## Déterminer quel échantillon analyser dans cette tâche

SAMPLE=\$(cat \$BASE_DIR/samplelist.txt | head -n \$SLURM_ARRAY_TASK_ID | tail -1)

## Créer un dossier pour cet échantillon et y naviguer
mkdir -p \$BASE_DIR/\$SAMPLE
cd \$BASE_DIR/\$SAMPLE

## Créer une liste avec le nom d'échantillon (nécessaire pour seqtk)
echo \"\$SAMPLE\" > tmp.list

## Boucle qui va chercher la séquence de l'échantillon de chaque locus et le mets dans un seul
## fichier .fasta, avec le nom de séquence qui inclut le nom du gène
for ALIGNMENT in $WD/$REF_TYPE/*$REF_FILE_ENDING
  do
    LOCUS=\$(basename \$ALIGNMENT $REF_FILE_ENDING)
	seqtk subseq \$ALIGNMENT tmp.list >> \$SAMPLE.fasta
	sed -i \"s/\$SAMPLE multi_.*/\$LOCUS/g\" \$SAMPLE.fasta
    sed -i \"s/\$SAMPLE single_.*/\$LOCUS/g\" \$SAMPLE.fasta
  done

## Supprimer la liste avec le nom d'échantillon
rm tmp.list

## Mapper les données Illumina de l'échantillon sur les séquences références de cet échantillon
bwa index \$SAMPLE.fasta
bwa mem -t 1 \$SAMPLE.fasta \\
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
$SRC/gatk-4.6.0.0/gatk HaplotypeCaller \\
  -R \$SAMPLE.fasta \\
  -I \$SAMPLE.rg.bam \\
  -O \$SAMPLE.vcf

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
## Soumettre ces les tâches de mappage et identification des variants
```bash
## Déterminer combien d'échantillons à analyser
NFILES=$(wc -l < samplelist2.txt)

sbatch --mail-user=$EMAIL --array=1-$NFILES bwa-gatk.sbatch
```
Job à resoumettre car terminée prématurément par manque de temps
Crataegus-ashei
Crataegus baroussana
Crataegus brainerdii
Crataegus coccinea coccinea
Crataegus coccinea pringlei
Crataegus cuprina
Crataegus gracilior
Crataegus harbisonii
Crataegus populnea
Crataegus rosei
Crataegus seratissima
Hesperomeles-obtusifolia_SRR12879568


## Supprimer les gènes ne contanant que des caractères indéterminée et les renommer avec l'extension filtered

```bash
nano filtre_fasta.py

##copier le code ci-dessous

#!/bin/bash

# Répertoire contenant les fichiers .fasta
repertoire_fasta="/scratch/mvallee/TP_session/Assem_avec_para/remap"

# Parcours de tous les fichiers .fasta dans les sous-dossiers
for fichier in "$repertoire_fasta"/*/*.fasta; do
    fichier_sortie="${fichier%.fasta}_filtered.fasta"
    
    # Initialiser une variable pour stocker la séquence en cours
    sequence=""
    header=""

    while IFS= read -r ligne; do
        # Si la ligne est un en-tête de séquence
        if [[ $ligne == ">"* ]]; then
            # Vérifier si la séquence précédente doit être écrite
            if [[ -n $sequence && $sequence != *[^N]* ]]; then
                # Si la séquence contient uniquement des N, on ne l'écrit pas
                sequence=""
            else
                # Si la séquence contient autre chose que des N, on l'écrit
                echo "$header" >> "$fichier_sortie"
                echo "$sequence" >> "$fichier_sortie"
            fi
            # Réinitialiser la séquence et enregistrer le nouvel en-tête
            header="$ligne"
            sequence=""
        else
            # Ajouter la ligne à la séquence en cours
            sequence+="$ligne"
        fi
    done < "$fichier"

    # Vérifier la dernière séquence après la boucle
    if [[ -n $sequence && $sequence != *[^N]* ]]; then
        sequence=""
    else
        echo "$header" >> "$fichier_sortie"
        echo "$sequence" >> "$fichier_sortie"
    fi
    
    echo "Fichier filtré créé : $fichier_sortie"
done
```

## Supprimer les N dans les séquences
```bash

# Répertoire contenant les fichiers _filtered.fasta
repertoire_fasta="/scratch/mvallee/TP_session/Assem_avec_para/remap"

# Parcours de tous les fichiers _filtered.fasta dans les sous-dossiers
for fichier in "$repertoire_fasta"/*/*_filtered.fasta; do
    # Créer un nouveau fichier de sortie en remplaçant le suffixe _filtered par _cleaned
    fichier_sortie="${fichier%.fasta}_cleaned.fasta"
    
    # Traite chaque ligne du fichier
    while IFS= read -r ligne; do
        # Si la ligne commence par '>', c'est une en-tête, on l'écrit telle quelle
        if [[ "$ligne" == \>* ]]; then
            echo "$ligne" >> "$fichier_sortie"
        else
            # Supprime les 'N' dans la séquence et l'écrit dans le fichier de sortie
            echo "${ligne//N/}" >> "$fichier_sortie"
        fi
    done < "$fichier"
    
    echo "Fichier nettoyé créé : $fichier_sortie"
done
```

# Créer des séquences avec hétérozygotes par gène
```bash
WD=/scratch/mvallee/TP_session/Assem_avec_para


## Créer dossier pour les listes de séquences par gène
mkdir -p $WD/remap/seqs_iupac
cd $WD/remap/seqs_iupac

## Créer une liste des gènes
ls /scratch/mvallee/TP_session/Assem_avec_para/exon_clean/*.FNA | xargs -n 1 basename | sed 's/\.FNA$//' > genelist.txt

## Pour chaque gène dans la liste de gène, aller chercher les séquences de chaque espèce et les combiner
## dans un seul fichier fasta
echo "while IGS= read -r GENE
  do
    echo \"Fetching sequences for \$GENE.\"
  while IGS= read -r SAMPLE
    do
      echo \"\$GENE\" > tmp.list
      seqtk subseq ../\$SAMPLE/\$SAMPLE.hetero.fasta tmp.list >> \$GENE.fasta
      sed -i \"s/\$GENE/\$SAMPLE/g\" \$GENE.fasta
    done < ../samplelist.txt
  done < genelist.txt" > fetch_seqs.sh
```
## Exécuter le shell
```bash
chmod +x fetch_seqs.sh
nohup ./fetch_seqs.sh > fetch_seqs.log &
```

## retirer les loci avec moins de 4 échantillons
```bash
WD=/scratch/mvallee/TP_session/Assem_avec_para

## Trouver les loci avec <4 échantillons et les supprimer
cd $WD/remap/seqs_iupac
for i in $(ls ./*.fasta)
  do
    NSAMPLES=$(grep '>' $i | wc -l)
    if [ $NSAMPLES -lt 4 ]
      then
        echo "$i contains $NSAMPLES samples, which is <4, so it was removed!"
        rm $i
    else
    echo "$i contains $NSAMPLES samples."
    fi
  done >> filter_min4.log
```
# Alignements avec MAFT
```bash
WD=/scratch/mvallee/TP_session/Assem_avec_para
EMAIL=marcoaurelevallee@gmail.com
TIME="0-12:00:00"

SEQS_PATH=$WD/remap/seqs_iupac/
TREES_PATH=$WD/remap/seqs_iupac/trees
mkdir -p $SEQS_PATH/align
mkdir -p $TREES_PATH/fastTree
cd $SEQS_PATH/align

## create a SLURM batchfile to align with MAFFT and run FastTree on each locus in an array
echo '#!/bin/bash' > mafft-fasttree.sbatch
echo "#SBATCH --job-name=mafft-fasttree
#SBATCH --output=mafft-fasttree-%a.out
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=$TIME

SEQIN=\$(ls $SEQS_PATH/*.fasta | head -n \$SLURM_ARRAY_TASK_ID | tail -1)
GENENAME=\$(basename \$SEQIN .fasta)

mafft --genafpair --maxiterate 1000 \$SEQIN > \$GENENAME.fasta

FastTree -nt -gtr \$GENENAME.fasta > $TREES_PATH/fastTree/\$GENENAME.tre" >> mafft-fasttree.sbatch
```
## Envoyer les tâches d'alignement à SLURM
```bash
## Déterminer le nombre de loci à analyser dans des sous-tâches séparées de SLURM
NFILES=$(ls -1 $SEQS_PATH/*.fasta | wc -l)

conda activate base
sbatch --mail-user=$EMAIL --array=1-$NFILES mafft-fasttree.sbatch
```

## Concatener les séquences alignées et les rendres lisibles
```bash
## Ajuster les variables ci-dessous de façon appropriée
SRC=/opt
WD=/scratch/mvallee/TP_session/Assem_avec_para
EMAIL=marcoaurelevallee@gmail.com
TIME="0-3:00:00"

## Créer dossier pour séquences concaténées
ALIGN_PATH=$WD/remap/seqs_iupac/align
mkdir -p $ALIGN_PATH/concat
cd $ALIGN_PATH/concat

## Créer un alignement concaténé format phyml avec catfasta2phyml
sbatch \
  --job-name=catfasta2phyml \
  --output=catfasta2phyml.log \
  --mail-user=$EMAIL \
  --nodes=1 \
  --time=$TIME \
  --cpus-per-task=1 \
  --mem-per-cpu=8G \
  --wrap="$SRC/catfasta2phyml/catfasta2phyml.pl \
  --concatenate $ALIGN_PATH/*.fasta \
  --sequential \
  1> concat.phy \
  2> concat.partitions"

## Créer un alignement concaténé format fasta avec catfasta2phyml
sbatch \
  --job-name=catfasta2phyml \
  --output=catfasta2phyml.log \
  --mail-user=$EMAIL \
  --nodes=1 \
  --time=$TIME \
  --cpus-per-task=1 \
  --mem-per-cpu=8G \
  --wrap="$SRC/catfasta2phyml/catfasta2phyml.pl \
  --fasta \
  --concatenate $ALIGN_PATH/*.fasta \
  1> concat.fasta \
  2> /dev/null"
```




## reformatter le fichier de partitions
```bash
ALIGN_PATH=$WD/remap/seqs_iupac/align/concat
cd $ALIGN_PATH
sed -i 's/^/DNA, /g' concat.partitions

## remove file path and ending from partition names
sed -i 's/DNA, .*\//DNA, /g' concat.partitions
sed -i 's/.fasta//g' concat.partitions
```

Télécharger l'alignement concatené sur l'ordinateur local et inspection visuelle
Les paralogues sont clairement mieux séparés dans l'alignement que par l'approche hybpiper+custom reference
Il y a quand même plusieurs zones assez divergentes donc je vais filtrer les alignements avec TAPER et Trimal puis refaire un concatenation


# Filtrer les alignements d'exons

```bash
ALIGNS_PATH=$WD/remap/seqs_iupac/align

mkdir -p $ALIGNS_PATH/taper
cd $ALIGNS_PATH

## Créer un fichier contenant la liste des alignements à filtrer avec TAPER
echo "" > align.list
for a in $ALIGNS_PATH/*.fasta
  do
    BASENAME=$(basename $a .fasta)
  echo "$a" >> align.list
  echo "$ALIGNS_PATH/taper/$BASENAME.fasta" >> align.list
  done
```
## créer un batchfile pour exécuter TAPER sur cette liste
```bash
echo '#!/bin/bash' > TAPER.sbatch
echo "#SBATCH --job-name=TAPER
#SBATCH --output=TAPER.out
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem-per-cpu=1G
#SBATCH --time=$TIME

## lancer TAPER sur la liste des gènes dans align.list
julia -t $THREADS $SRC/TAPER-1.0.0/correction_multi.jl \\
  --mask $MASKING_CHAR \\
  --cutoff $FILTERING_CUTOFF \\
  --list align.list" >> TAPER.sbatch
```

## envoyer la tâche TAPER à SLURM
```bash
sbatch --mail-user=$EMAIL TAPER.sbatch
```

## Filtration avec Trimal

```bash
SRC=/opt
WD=/scratch/mvallee/TP_session/Assem_avec_para
EMAIL=votre.courriel@umontreal.ca
TIME="0-3:00:00"

## Ajuster aussi les paramètres de trimAl si vous voulez
WINDOW_SIZE=3
MIN_NONGAP_PERCENT=0.5
MIN_SIMILARITY=0.7
MIN_OVERLAP=0.5

ALIGNS_PATH=$WD/remap/seqs_iupac/align

mkdir -p $ALIGNS_PATH/trimal
cd $ALIGNS_PATH

## créer un batchfile pour filtrer les alignements
echo '#!/bin/bash' > trimal.sbatch
echo "#SBATCH --job-name=trimal
#SBATCH --output=trimal-%a.out
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=$TIME

## sélectionner l'alignement à filtrer dans cette tâche
ALIGNIN=\$(ls $ALIGNS_PATH/taper/*.fasta | head -n \$SLURM_ARRAY_TASK_ID | tail -1)

## déterminer le nom du gène
GENENAME=\$(basename \$ALIGNIN .fasta)

## lancer trimAl
$SRC/trimAl_1.5.0/trimal \\
  -in \$ALIGNIN \\
  -out ./trimal/\$GENENAME.fasta \\
  -w $WINDOW_SIZE \\
  -gapthreshold $MIN_NONGAP_PERCENT \\
  -simthreshold $MIN_SIMILARITY \\
  -resoverlap $MIN_OVERLAP \\
  -seqoverlap $MIN_OVERLAP" >> trimal.sbatch
```

## exécuter trimAl sur tous les alignements filtrés par TAPER
```bash
NFILES=$(ls -1 $ALIGNS_PATH/taper/*.fasta | wc -l)
sbatch --mail-user=$EMAIL --array=1-$NFILES trimal.sbatch
```

# Concatener les séquences avec code IUPAC
```bash
SRC=/opt
WD=/scratch/mvallee/TP_session/Assem_avec_para
EMAIL=marcoaurelevallee@gmail.com
TIME="0-3:00:00"

## Créer dossier pour séquences concaténées
ALIGN_PATH=$WD/remap/seqs_iupac/align/trimal
mkdir -p $ALIGN_PATH/concat
cd $ALIGN_PATH/concat

## Créer un alignement concaténé format phyml avec catfasta2phyml
sbatch \
  --job-name=catfasta2phyml \
  --output=catfasta2phyml.log \
  --mail-user=$EMAIL \
  --nodes=1 \
  --time=$TIME \
  --cpus-per-task=1 \
  --mem-per-cpu=8G \
  --wrap="$SRC/catfasta2phyml/catfasta2phyml.pl \
  --concatenate $ALIGN_PATH/*.fasta \
  --sequential \
  1> concat.phy \
  2> concat.partitions"

## Créer un alignement concaténé format fasta avec catfasta2phyml
sbatch \
  --job-name=catfasta2phyml \
  --output=catfasta2phyml.log \
  --mail-user=$EMAIL \
  --nodes=1 \
  --time=$TIME \
  --cpus-per-task=1 \
  --mem-per-cpu=8G \
  --wrap="$SRC/catfasta2phyml/catfasta2phyml.pl \
  --fasta \
  --concatenate $ALIGN_PATH/*.fasta \
  1> concat.fasta \
  2> /dev/null"

```

## Supprimer les colonnes contenant presque uniquement des caractères IUPAC

```bash
# Définir les caractères IUPAC
IUPAC="RYWSKMBDHVN"

# Fichier d'entrée et sortie
input_file="concat.phy"
output_file="concat-trimmed.phy"

# Lire le fichier et traiter ligne par ligne
awk '{print $0}' $input_file | while read line; do
    # Si ce n'est pas une ligne d'en-tête (la première ligne de l'alignement), on continue
    if [[ ! "$line" =~ ^[0-9]+$ ]]; then
        seq=$line
        seq_length=${#seq}
        # Créer une nouvelle chaîne pour la séquence filtrée
        new_seq=""
        
        # Traiter chaque colonne de la séquence
        for (( i=0; i<$seq_length; i++ )); do
            char="${seq:i:1}"
            # Vérifier si le caractère est un IUPAC
            if [[ $IUPAC == *"$char"* ]]; then
                new_seq+="-"
            else
                new_seq+="$char"
            fi
        done
        
        # Ajouter la séquence modifiée au fichier de sortie
        echo $new_seq >> $output_file
    fi
done
```

## Reformater le fichier de partitions
```bash
ALIGN_PATH=$WD/remap/seqs_iupac/align/concat
cd $ALIGN_PATH

## add the "DNA, " required at start of each line of partition file for RAxML
sed -i 's/^/DNA, /g' concat.partitions

## remove file path and ending from partition names
sed -i 's/DNA, .*\//DNA, /g' concat.partitions
sed -i 's/.fasta//g' concat.partitions
```

Le flitrage a bien fonctionnée, c'Est beaucoup mieux!!

# Phasage des allèles par polarisation itérative avec IterPol
```bash
SRC=/opt
WD=/scratch/mvallee/TP_session/Assem_avec_para/remap/seqs_iupac
ALIGNMENT=/scratch/mvallee/TP_session/Assem_avec_para/remap/seqs_iupac/align/trimal/concat/concat_trim.phy
mkdir -p $WD/iterpol
cd $WD/iterpol

sbatch --job-name=IterPol \
  --output=IterPol.out \
  --mail-user=marcoaurelevallee@gmail.com \
  --mail-type=end \
  --time=48:00:00 \
  --wrap "source activate bio; python ~/TP_session/ref/IterPol_v0.4.py \
    --infile $ALIGNMENT \
    --out_prefix crataegus \
    --method raxml \
    --threads 20 \
    --raxml_path \"$SRC/RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3\" > output.log 2>&1"
```

## Analyse de Maximum de vraisemblance dans RAxML avec alignements polarisés
```bash
SRC_RAXML=/opt/RAxML-8.2.12
WD=/scratch/mvallee/TP_session/Assem_avec_para/remap/seqs_iupac/iterpol/raxml
ALIGNMENT=/scratch/mvallee/TP_session/Assem_avec_para/remap/seqs_iupac/iterpol/Crataegus_final.polarized.phy
MODEL=GTRCAT
OUTPUT_PREFIX=Crataegus
BOOTSTRAP_REPS=1000
EMAIL=marcoaurelevallee@gmail.com
TIME="0-50:00:00"


cd $WD

ln $ALIGNMENT $OUTPUT_PREFIX.phy

## Exécuter ce batchfile en mode non-intéractif sur SLURM
sbatch \
  --job-name=raxml \
  --output=raxml.log \
  --mail-user=$EMAIL \
  --nodes=1 \
  --time=$TIME \
  --cpus-per-task=6 \
  --mem-per-cpu=2G \
  --wrap="$SRC_RAXML/raxmlHPC-PTHREADS-SSE3 \
    $ARGUMENTS_FOR_RAXML \
    -f a \
    -T 6 \
    -p 1234 \
    -x 1234 \
    -m $MODEL \
    -# $BOOTSTRAP_REPS \
    -s Crataegus.phy.reduced \
    -n $OUTPUT_PREFIX"

```


# Arbre d'espèces avec ASTRAL
```bash
## Générer des arbres de gènes avec les données filtré 
WD=/scratch/mvallee/TP_session/Assem_avec_para
EMAIL=marcoaurelevallee@gmail.com
TIME="0-12:00:00"

SEQS_PATH=/scratch/mvallee/TP_session/Assem_avec_para/remap/seqs_iupac/align/trimal
TREES_PATH=$WD/remap/seqs_iupac/trees
cd $TREES_PATH/fastTree
echo '#!/bin/bash' > fasttree-trim.sbatch
echo "#SBATCH --job-name=fasttree-trim
#SBATCH --output=fasttree-trim-%a.out
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=$TIME

SEQIN=\$(ls $SEQS_PATH/*.fasta | head -n \$SLURM_ARRAY_TASK_ID | tail -1)
GENENAME=\$(basename \$SEQIN .fasta)


FastTree -nt -gtr \$GENENAME.fasta > $TREES_PATH/fastTree/\$GENENAME.tre" >> fasttree-trim.sbatch
```
## Exécuter FastTree
```bash
## Déterminer le nombre de loci à analyser dans des sous-tâches séparées de SLURM
NFILES=$(ls -1 $SEQS_PATH/*.fasta | wc -l)

## Envoyer les tâches d'alignement à SLURM
conda activate base
sbatch --mail-user=$EMAIL --array=1-$NFILES fasttree-trim.sbatch
```


tous les fichiers .tre du répertoire /scratch/mvallee/TP_session/Assem_avec_para/remap/seqs_iupac/trees/fastTree 
ont été transféré sur l'ordinateur local puis dans le cluster beluga

## Arbre d'espèces avec ASTRAL
Sur le cluster Beluga
```bash
SRC_ASTRAL=/project/def-bourret/shared/progs/Astral
WD=/home/maval/projects/def-bourret/maval/Seqs_iupac/trees/Astral
GENETREES=/home/maval/projects/def-bourret/maval/Seqs_iupac/trees/fastTree-trim/*.tre
OUTGROUP=Malus-domestica_Mlu855Rosi
EMAIL=marcoaurelevallee@gmail.com
THREADS=16

mkdir -p $WD
cd $WD

## create a list of all gene trees we want to input to astral
cat $GENETREES > gene.trees

## check how many gene trees we have in this file
wc -l gene.trees

## run astral
module load StdEnv/2020
module load java/17.0.2

sbatch --job-name=Astral \
  --output=Astral.out \
  --mail-user=$EMAIL \
  --mail-type=end \
  --time=04:00:00 \
  --mem-per-cpu=1G \
  --cpus-per-task=$THREADS \
  --wrap "java -jar $SRC_ASTRAL/astral.*.jar --input gene.trees --output astral.tre --outgroup $OUTGROUP --branch-annotate 2"
```


## Créer un fichier .tre avec seulement les valeurs pp1
```bash
input_file="astral.tre"  
output_file="astral_pp1.tre"

# Lire chaque ligne du fichier d'entrée
while IFS= read -r line; do
    
 modified_line=$(echo "$line" | sed -E 's/\[([^]]*pp1=([^;]*))[^]]*\]/[pp1=\2]/g')
    # Ajouter la ligne modifiée au fichier de sortie
    echo "$modified_line" >> "$output_file"
done < "$input_file"

echo "Extraction des valeurs pp1 terminée. Résultat enregistré dans $output_file."
```
