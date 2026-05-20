
# Workshop Analyse Bioinformatique sous Galaxy 2026 — 

---

# 📚 Sommaire

* [Introduction](#-introduction)
* [Workflow bioinformatique](#-workflow-bioinformatique)
* [QC des données](#qc-des-donnees)
* [Assemblage génomique](#assemblage-genomique)
* [Recherche de variants](#recherche-de-variants)
* [Analyse des AMRs](#analyse-des-amrs)
* [Notion de Workflow](#notion-de-workflow)
* [Tutoriel sur sars-cov2](#tutoriel-sur-sars-cov2)

---

# 🔬 Introduction

Bienvenue dans cet atelier d’analyse bioinformatique.

Ce support permet aux participants de :

* suivre les étapes d’un workflow bioinformatique ;
* Copier les liens d'acces aux données;
* cliquer pour afficher les réponses ;
* apprendre progressivement les concepts clés.

---

# 🧪 Workflow bioinformatique

```text
FASTQ brut
    ↓
Contrôle qualité (FastQC)
    ↓
Trimming / nettoyage
    ↓
Alignement / assemblage Recherche de variants / Annotation Recherche de Gene AMR
    ↓
Analyse biologique
```

---


# QC des donnees
> Basé sur le tutoriel GTN : [gxy.io/GTN:T00239](https://gxy.io/GTN:T00239

Lors du séquençage, les bases nucléotidiques d'un échantillon d'ADN ou d'ARN (bibliothèque) sont déterminées par le séquenceur. Pour chaque fragment de la bibliothèque, une séquence est générée, également appelée lecture , qui est simplement une succession de nucléotides.

Les technologies de séquençage modernes permettent de générer un grand nombre de séquences en une seule expérience. Cependant, aucune technologie n'est parfaite et chaque instrument produit des erreurs de nature et de quantité variables, comme l'identification incorrecte de nucléotides. Ces erreurs d'identification sont dues aux limitations techniques propres à chaque plateforme de séquençage.

Il est donc nécessaire de comprendre, d'identifier et d'éliminer les types d'erreurs susceptibles d'affecter l'interprétation des analyses ultérieures. Le contrôle qualité des séquences constitue ainsi une première étape essentielle de votre analyse. Détecter les erreurs au plus tôt permet de gagner du temps par la suite.

## Téléchargé un fichier de séquence brute
1. créez un nouveau historique
2. renommez le
3. importer les données NGS via Faster Download and Extract Reads in FASTQ en utilisant le ID SRR suivant 
````
SRR3111247
````

## Lancer l'outil FATSQC sur ces données
1. commbien des reads nous avons sur ces données?
2. Quelle la taille des reads?
3. Qu'est ce qu'il faut faire pour améliorer  ce jeux de données?
## Lancer fastp sur ce jeux de données en utlisant un cut-off de 50 pb en taille.
1. comparer le resumltats de filtrage avec le resultat de FASTQ.

---
# 🧬 Assemblage génomique

 



---

## Recherche de variants



---
## Analyse des AMRs


---
## Analyse des AMRs


---
## Notion de Workflow


---
## Tutoriel sur sars-cov2



---

