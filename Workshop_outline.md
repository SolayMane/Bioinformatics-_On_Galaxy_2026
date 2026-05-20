
# Workshop Analyse Bioinformatique sous Galaxy 2026 — 

---

# 📚 Sommaire

* [Introduction](#-introduction)
* [Workflow bioinformatique](#-workflow-bioinformatique)
* [QC des données](#qc-des-donnees)
* [Assemblage génomique](#-assemblage-genomique)
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
Alignement / assemblage Recherche de variants / Annotation / Recherche de Gene AMR
    ↓
Analyse biologique
```

---


# QC des donnees


Lors du séquençage, les bases nucléotidiques d'un échantillon d'ADN ou d'ARN (bibliothèque) sont déterminées par le séquenceur. Pour chaque fragment de la bibliothèque, une séquence est générée, également appelée lecture , qui est simplement une succession de nucléotides.

Les technologies de séquençage modernes permettent de générer un grand nombre de séquences en une seule expérience. Cependant, aucune technologie n'est parfaite et chaque instrument produit des erreurs de nature et de quantité variables, comme l'identification incorrecte de nucléotides. Ces erreurs d'identification sont dues aux limitations techniques propres à chaque plateforme de séquençage.

Il est donc nécessaire de comprendre, d'identifier et d'éliminer les types d'erreurs susceptibles d'affecter l'interprétation des analyses ultérieures. Le contrôle qualité des séquences constitue ainsi une première étape essentielle de votre analyse. Détecter les erreurs au plus tôt permet de gagner du temps par la suite.

## Téléchargé un fichier de séquence brute à partir de la base de données SRA
1. créez un nouveau historique
2. renommez le
3. importer les données NGS via Faster Download and Extract Reads in FASTQ en utilisant le ID SRR suivant 
````
SRR3111247
````
4. inspecter un des ces fichier
````
@SRR3111247.1/1
GGAATGCCTGATGGCGGTTCGGCACCTGGTTTGCTGAGAGACATCGCTCGCTGCGCATACCATGACGAATAGGGACTGTCGCGGTATGCGTTGCTGCTAA
+
B>:A<9>A@9(BBB95@BD@?C@C????ACCA>>CCDC@>9@A:>>CCFEGHJIHEA>CFC9GGHIJJIGGGD9?GHEIBFIJIIIIHHHDDDFFFFCC@
````

Qu'est ce que signifie le score de qualité ?

Le score de qualité de chaque séquence est une chaîne de caractères, un pour chaque base de la séquence nucléotidique, utilisée pour caractériser la probabilité d'une identification erronée de chaque base. Ce score est encodé à l'aide de la table de caractères ASCII (avec quelques différences historiques ).

Pour économiser de l'espace, le séquenceur enregistre un caractère ASCII pour représenter les scores de 0 à 42. Par exemple, 10 correspond à « + » et 40 à « I ». FastQC sait interpréter ces données. On parle alors de score « Phred ».

Encodage du score de qualité en caractères ASCII pour différents encodages Phred. La séquence de code ASCII est affichée en haut, avec les symboles de 33 à 64, les lettres majuscules, d'autres symboles, puis les lettres minuscules. Sanger couvre la plage de 33 à 73, tandis que Solexa est décalé, de 59 à 104. Illumina 1.3 commence à 54 et va jusqu'à 104. Illumina 1.5 est décalé de trois scores vers la droite, mais se termine toujours à 104. Illumina 1.8+ reprend la plage de Sanger, avec un score de plus.

 <img width="3436" height="656" alt="image" src="https://github.com/user-attachments/assets/7cb40a95-79f0-4ed8-9106-6ed4470067b9" />

Chaque nucléotide est donc associé à un caractère ASCII représentant son score de qualité Phred , c'est-à-dire la probabilité d'une identification de base incorrecte :


| Score de qualité Phred | Probabilité d'appel de base incorrect | Précision de l'appel de base |
| :---: | :--- | :--- |
| 10 | 1 sur 10 | 90% |
| 20 | 1 sur 100 | 99% |
| 30 | 1 sur 1000 | 99,9% |
| 40 | 1 sur 10 000 | 99,99% |
| 50 | 1 sur 100 000 | 99,999% |
| 60 | 1 sur 1 000 000 | 99,9999% |

Que représentent les valeurs de 0 à 42 ? Ces nombres, intégrés dans une formule, nous indiquent la probabilité d’une erreur pour cette base. Voici la formule, où Q est notre score de qualité (0 à 42) et P la probabilité d’une erreur :
```
Q = -10 log10(P)
```
En utilisant cette formule, nous pouvons calculer qu'un score de qualité de 40 signifie seulement 0,00010 de probabilité d'erreur !

Question: Comment calculer la précision de la base nucléotidique avec le code ASCII ```/```?

Cela peut être calculé comme suit

**Le code ASCII /est 47**

**score de qualité = 47-33=14**

**Formule pour calculer la probabilité d'erreur :P =10− Q / 10**

**Probabilité d'erreur =10− 14 / 10 = 0,03981 = 3,981 %**

**Par conséquent, la précision est de 100 - 3,981 = 96,019 %.**


## Lancer l'outil FATSQC sur ces données
1. commbien des reads nous avons sur ces données?
2. Quelle la taille des reads?
3. Qu'est ce qu'il faut faire pour améliorer  ce jeux de données?
## Lancer fastp sur ce jeux de données en utlisant ces paramètres :
* cut-off de 50 pb en taille.
* Phred score de 20
1. comparer le resultat de filtrage avec le resultat précédant de FASTQC.

## Téléchargé des données à partir d'un lien
1. importer les données NGS en utilsant ces liens
````
https://zenodo.org/record/3977236/files/female_oral2.fastq-4143.gz
````
2. Lancer FASTQC et inspecter le resultats
3. Comment vous trouvez ces données en terme de qualité
---
# 🧬Assemblage génomique
Dans ce tuotriel nous allons faire un assemblage d'un seul genome puis nous allons l'annoter puis le visualiser
 1. Creer un nouveau historique
 2. Télécharger les données brutes ces leins :
````
https://zenodo.org/record/10669812/files/DRR187559_1.fastqsanger.bz2
https://zenodo.org/record/10669812/files/DRR187559_2.fastqsanger.bz2
````



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

