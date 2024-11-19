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
```
