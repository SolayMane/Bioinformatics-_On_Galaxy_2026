# 🧬 Workshop Bioinformatique sous Galaxy — 2026

---

# 📚 Sommaire

- [Introduction générale](#-partie-i--introduction-générale)
- [Contrôle qualité des données NGS](#-partie-ii--contrôle-qualité-des-données-ngs)
- [Assemblage génomique bactérien](#-partie-iii--assemblage-génomique-bactérien)
- [Annotation génomique bactérienne](#-partie-iv--annotation-génomique-bactérienne)
- [Recherche de variants](#-partie-v--recherche-de-variants)
- [Détection des gènes AMR](#-partie-vi--détection-des-gènes-amr)
- [Analyse SARS-CoV-2](#-partie-vii--analyse-sars-cov-2)
- [Notions avancées](#-partie-viii--notions-avancées)

---

# 🧬 Partie I — Introduction générale

## 🔬 Introduction à la bioinformatique

La bioinformatique est une discipline qui combine :

- la biologie ;
- l’informatique ;
- les statistiques ;
- les mathématiques.

Elle permet d’analyser les données biologiques générées par les technologies modernes de séquençage.

### Applications principales

- Génomique
- Transcriptomique
- Métagénomique
- Épidémiologie génomique
- Détection des résistances aux antibiotiques
- Analyse comparative de génomes

---

## 🌌 Introduction à Galaxy

Galaxy est une plateforme bioinformatique accessible via navigateur web permettant :

- d’exécuter des analyses bioinformatiques ;
- de construire des workflows reproductibles ;
- de partager des historiques ;
- d’éviter l’utilisation intensive de la ligne de commande.

### Avantages de Galaxy

- Interface graphique intuitive
- Reproductibilité
- Historique complet des analyses
- Partage des workflows
- Compatible avec de nombreux outils bioinformatiques

---

## 🧪 Workflow bioinformatique général

```text
                ┌───────────────────┐
                │   Données FASTQ   │
                └─────────┬─────────┘
                          │
                          ▼
                ┌───────────────────┐
                │ Contrôle Qualité  │
                │   FastQC / Falco  │
                └─────────┬─────────┘
                          │
                          ▼
                ┌───────────────────┐
                │ Nettoyage Reads   │
                │ fastp / cutadapt  │
                └─────────┬─────────┘
                          │
            ┌─────────────┴─────────────┐
            ▼                           ▼
┌───────────────────┐      ┌────────────────────┐
│ Alignement Reads  │      │ Assemblage Genome  │
│ bwa / bowtie2     │      │ Shovill / SPAdes   │
└─────────┬─────────┘      └─────────┬──────────┘
          │                           │
          ▼                           ▼
┌───────────────────┐      ┌────────────────────┐
│ Variant Calling   │      │ Annotation Genome  │
│ Snippy / Freebayes│      │ Prokka / Bakta     │
└─────────┬─────────┘      └─────────┬──────────┘
          │                           │
          └─────────────┬─────────────┘
                        ▼
             ┌────────────────────┐
             │ Visualisation      │
             │ JBrowse / Circos   │
             └─────────┬──────────┘
                       ▼
             ┌────────────────────┐
             │ Interprétation     │
             │ biologique         │
             └────────────────────┘
```

---

# 🧪 Partie II — Contrôle qualité des données NGS

## 📖 Introduction aux fichiers FASTQ

Lors du séquençage, les bases nucléotidiques d'un échantillon d'ADN ou d'ARN sont déterminées par le séquenceur. Chaque fragment génère une lecture appelée **read**.

Les technologies modernes permettent de produire des millions de reads, mais aucune plateforme n’est parfaite.

Le contrôle qualité constitue donc une étape essentielle avant toute analyse.

---

## 📥 Téléchargement des données depuis SRA

### Base de données SRA

La base SRA (Sequence Read Archive) est le plus grand dépôt mondial de données NGS.

Elle est synchronisée entre :

- NCBI SRA
- ENA
- DRA

---

## 🚀 Étapes sous Galaxy

1. Créer un nouvel historique
2. Renommer l’historique
3. Explorer SRA
4. Utiliser :

```text
Faster Download and Extract Reads in FASTQ
```

### Accession utilisée

```text
SRR3111247
```

---

## 📄 Structure d’un fichier FASTQ

```text
@SRR3111247.1/1
GGAATGCCTGATGGCGGTTCGGCACCTGGTTTGCTGAGAGACATCGCTCGCTGCGCATACCATGACGAATAGGGACTGTCGCGGTATGCGTTGCTGCTAA
+
B>:A<9>A@9(BBB95@BD@?C@C????ACCA>>CCDC@>9@A:>>CCFEGHJIHEA>CFC9GGHIJJIGGGD9?GHEIBFIJIIIIHHHDDDFFFFCC@
```

---

## 📊 Scores Phred et qualité des reads

Le score Phred représente la probabilité d’erreur d’identification d’une base.

### Formule

```text
Q = -10 log10(P)
```

### Tableau de qualité

| Score Phred | Probabilité d'erreur | Précision |
|---|---|---|
| 10 | 1/10 | 90% |
| 20 | 1/100 | 99% |
| 30 | 1/1000 | 99.9% |
| 40 | 1/10000 | 99.99% |

---

<img width="3436" height="656" alt="image" src="https://github.com/user-attachments/assets/7cb40a95-79f0-4ed8-9106-6ed4470067b9" />

---

## 🧪 Contrôle qualité avec FastQC

### Objectifs

Détecter :

- mauvaise qualité ;
- adaptateurs ;
- biais GC ;
- contamination ;
- reads trop courts.

---

## 🚀 Lancer FastQC

<img width="293" height="458" alt="image" src="https://github.com/user-attachments/assets/66599f75-5cfc-4607-83fb-57bf3c2ec0e3" />

---

## ❓ Questions

1. Quel encodage Phred est utilisé ?
2. Combien de reads sont présents ?
3. Quelle est la taille moyenne des reads ?

---

## ✂️ Nettoyage avec fastp

### Paramètres

- Taille minimale : 50 pb
- Score Phred minimal : 20

<img width="297" height="279" alt="image" src="https://github.com/user-attachments/assets/0143828a-e481-4307-a6df-e5c29516a3d8" />

---

## 🔄 Workflow QC complet

```text
FASTQ brut
     │
     ▼
FastQC / Falco
     │
     ▼
Détection :
- mauvaise qualité
- adaptateurs
- biais GC
     │
     ▼
fastp / cutadapt
     │
     ▼
FASTQ nettoyé
     │
     ▼
FastQC final
     │
     ▼
Validation qualité
```

---

# 🧬 Partie III — Assemblage génomique bactérien

# 🦠 Cas d'étude 1 — Assemblage d’un génome MRSA

---

## 🎯 Objectifs

À la fin de ce module :

- importer des reads paired-end ;
- nettoyer les données ;
- assembler un génome ;
- évaluer l’assemblage ;
- visualiser les résultats.

---

## 🧬 Contexte biologique

Les données proviennent de :

```text
Methicillin-Resistant Staphylococcus aureus (MRSA)
```

Applications :

- recherche de gènes AMR ;
- épidémiologie ;
- comparaison de souches ;
- typage génomique.

---

## 📥 Importation des données

```text
https://zenodo.org/record/10669812/files/DRR187559_1.fastqsanger.bz2
https://zenodo.org/record/10669812/files/DRR187559_2.fastqsanger.bz2
```

---

## 🧪 Contrôle qualité des reads

Outil :

```text
Falco ou FastQC
```

Analyser :

- Per base quality
- GC content
- Adapter content

---

## ✂️ Nettoyage avec fastp

### Paramètres

```text
Length required: 30
Cutting mean quality: 20
```

---

## 🧩 Assemblage avec Shovill

### Pourquoi Shovill ?

Shovill est optimisé pour les petits génomes bactériens.

Basé sur :

```text
SPAdes
```

---

## 📤 Résultats générés

- contigs FASTA ;
- log ;
- graphe d’assemblage.

---

## 📊 Évaluation avec QUAST

Métriques importantes :

- N50 ;
- nombre de contigs ;
- GC ;
- taille totale.

---

## 🧬 Visualisation avec JBrowse

<img width="298" height="392" alt="image" src="https://github.com/user-attachments/assets/33e045de-e780-4e19-b956-849623ca6bb7" />

<img width="1303" height="445" alt="image" src="https://github.com/user-attachments/assets/fb41909c-7b33-4edc-bbe6-f1c0dcd9020e" />

---

## 🌌 Visualisation avec Circos

### Objectif

Créer une représentation circulaire :

- des gènes ;
- de leur densité ;
- des régions génomiques.

---

## 🔄 Workflow Assemblage complet

```text
Reads Illumina
       │
       ▼
FastQC / Falco
       │
       ▼
fastp
       │
       ▼
Shovill
       │
       ▼
Contigs
       │
       ▼
QUAST
       │
       ├──────────────┐
       ▼              ▼
JBrowse            Circos
       │              │
       └──────┬───────┘
              ▼
     Interprétation biologique
```

---

# 🦠 Cas d'étude 2 — Assemblage multi-génomes Salmonella

## 🧬 Contexte biologique

Les données proviennent de :

```text
Salmonella gallinarum
```

Applications :

- comparaison de souches ;
- étude épidémiologique ;
- AMR ;
- annotation comparative.

---

## 📥 Importation de l’historique partagé

```text
https://usegalaxy.eu/u/solay/h/inh-workshop-2026
```

<img width="905" height="500" alt="image" src="https://github.com/user-attachments/assets/3485978e-26ae-4bb5-b253-33c32a958cd2" />

---

## 🧪 Étapes du workflow

1. QC
2. Nettoyage
3. Assemblage
4. Annotation
5. Comparaison génomique

---

# 🧬 Partie IV — Annotation génomique bactérienne

## 🧠 Introduction à l’annotation

L’annotation génomique consiste à identifier :

- les gènes ;
- les ARNr ;
- les ARNt ;
- les protéines.

---

## 🚀 Annotation avec Prokka

### Outil

```text
Prokka
```

---

## 📥 Entrée

```text
Contigs assemblés
```

---

## 📂 Fichiers générés

### `.gff`

Fichier principal contenant :

- coordonnées ;
- gènes ;
- annotations.

### `.faa`

Séquences protéiques.

### `.ffn`

Séquences nucléotidiques des gènes.

### `.gbk`

Format GenBank.

### `.sqn`

Soumission NCBI.

---

## 🔄 Workflow Annotation

```text
Contigs assemblés
        │
        ▼
      Prokka
        │
 ┌──────┼────────┐
 ▼      ▼        ▼
GFF    FAA      GBK
 │       │        │
 ▼       ▼        ▼
JBrowse BLAST  Soumission NCBI
```

---

# 🧬 Partie V — Recherche de variants

## 🧬 Introduction

Le variant calling permet d’identifier :

- SNPs ;
- insertions ;
- délétions.

---

## 🧪 Étapes principales

1. Alignement
2. Tri BAM
3. Indexation
4. Variant calling
5. Génération VCF

---

## 🛠️ Outils possibles

- bwa
- bowtie2
- samtools
- Snippy
- Freebayes

---

## 🔄 Workflow SNP

```text
Reads
   │
   ▼
Alignement référence
   │
   ▼
BAM trié/indexé
   │
   ▼
Variant Calling
   │
   ▼
VCF
   │
   ├─────────────┐
   ▼             ▼
Annotation     Visualisation
SNP            IGV/JBrowse
```

---

# 🦠 Partie VI — Détection des gènes AMR

## 🔬 Introduction

La résistance aux antibiotiques est un problème majeur de santé publique.

---

## 🎯 Objectifs

- détecter les ARGs ;
- identifier les plasmides ;
- visualiser les gènes AMR.

---

## 🛠️ Outil principal

```text
staramr
```

---

## 📂 Résultats générés

| Fichier | Description |
|---|---|
| summary.tsv | Résumé |
| detailed_summary.tsv | Détails |
| plasmidfinder.tsv | Plasmides |

---

## ❓ Gène attendu chez MRSA

```text
mecA
```

---

## 🧬 Conversion en GFF3

Outil :

```text
Table to GFF3
```

---

## 🧬 Annotation avec Bakta

Bakta permet :

- annotation rapide ;
- standardisation ;
- visualisation génomique.

---

## 🌌 Visualisation avec JBrowse

Visualiser :

- ARGs ;
- plasmides ;
- annotations.

---

## 🔄 Workflow AMR

```text
Contigs assemblés
        │
        ▼
      staramr
        │
        ├─────────────┐
        ▼             ▼
 ARGs détectés     Plasmides
        │
        ▼
 Conversion GFF3
        │
        ▼
      Bakta
        │
        ▼
     JBrowse
        │
        ▼
Interprétation AMR
```

---

# 🦠 Partie VII — Analyse SARS‑CoV‑2

## 🎯 Objectifs

- alignement viral ;
- variant calling ;
- génération de consensus ;
- assignation de lignages.

---

## 🛠️ Outils possibles

- bwa
- ivar
- pangolin
- nextclade

---

## 🔄 Workflow SARS‑CoV‑2

```text
Reads viraux
      │
      ▼
QC + trimming
      │
      ▼
Alignement Wuhan
      │
      ▼
Variant Calling
      │
      ▼
Consensus FASTA
      │
      ▼
Pangolin / Nextclade
      │
      ▼
Interprétation épidémiologique
```

---

# 📚 Partie VIII — Notions avancées

## 🧠 Reproductibilité

Galaxy permet :

- sauvegarde d’historiques ;
- partage ;
- workflows automatiques.

---

## 🔄 Notion de Workflow

Un workflow permet :

- automatisation ;
- standardisation ;
- reproductibilité.

---

## 📦 Export des workflows

Formats :

```text
.ga
```

---

# 📁 Structure recommandée du dépôt GitHub

```text
galaxy-bioinformatics-workshop/
│
├── README.md
├── docs/
│   ├── images/
│   ├── workflows/
│   └── tutorials/
│
├── data/
├── slides/
├── workflows-ga/
└── environment/
```

---

# ✅ Conclusion

Ce workshop introduit :

- l’analyse NGS sous Galaxy ;
- le contrôle qualité ;
- l’assemblage bactérien ;
- l’annotation génomique ;
- la détection AMR ;
- le variant calling ;
- la visualisation génomique.

L’objectif est de fournir une approche pratique, reproductible et pédagogique de la bioinformatique moderne.

