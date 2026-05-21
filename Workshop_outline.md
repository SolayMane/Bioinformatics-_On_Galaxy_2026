
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
En bioinformatique, la base de données SRA (Sequence Read Archive) est le plus grand référentiel public mondial de données de séquençage de l'ADN et de l'ARN à haut débit.
Elle fait partie intégrante de l'INSDC (International Nucleotide Sequence Database Collaboration). Les données soumises à l'une des plateformes suivantes sont automatiquement synchronisées entre elles :
NCBI SRA (États-Unis)
ERA / ENA (European Nucleotide Archive, gérée par l'EBI en Europe)
DRA (DDBJ Sequence Read Archive, au Japon)

1. Créez un nouveau historique
2. Renommez le
3. Comment explorer la base de données SRA (Demo NCBI)
4. Importer les données NGS via ````Faster Download and Extract Reads in FASTQ```` en utilisant 'l'accession number' SRRXXX suivant 
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
<detailS>
<summary>Qu'est ce que signifie le score de qualité ?</summary>

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

Que représentent les valeurs de 0 à 42 ?
Ces nombres, intégrés dans une formule, nous indiquent la probabilité d’une erreur pour cette base. Voici la formule, où Q est notre score de qualité (0 à 42) et P la probabilité d’une erreur :
```
Q = -10 log10(P)
```
En utilisant cette formule, nous pouvons calculer qu'un score de qualité de 40 signifie seulement 0,00010 de probabilité d'erreur !

Question: Comment calculer la précision de la base nucléotidique avec le code ASCII ```/```?

Cela peut être calculé comme suit

Le code ASCII /est 47

score de qualité = 47-33=14

Formule pour calculer la probabilité d'erreur :P =10−<sup>Q/10</sup>

Probabilité d'erreur =10− 14 / 10 = 0,03981 = 3,981 %

Par conséquent, la précision est de 100 - 3,981 = 96,019 %.

</details>

### Évaluez la qualité avec FastQC
Une autre méthode pour vérifier la qualité des séquences consiste à utiliser FastQC . Cet outil propose un ensemble d'analyses modulaires permettant de détecter d'éventuels problèmes dans vos données avant toute analyse ultérieure. 

1.Lancer ````FastQC````

<img width="293" height="458" alt="image" src="https://github.com/user-attachments/assets/66599f75-5cfc-4607-83fb-57bf3c2ec0e3" />


2. Quel encodage Phred est utilisé dans le fichier FASTQ pour ces séquences ?
3. commbien des reads nous avons sur ces données?
4. Quelle la taille des reads?
5. Lancer `fastp` sur ce jeux de données en utlisant ces paramètres :
* cut-off de 50 pb en taille.
* Phred score de 20
  
<img width="297" height="279" alt="image" src="https://github.com/user-attachments/assets/0143828a-e481-4307-a6df-e5c29516a3d8" />

6. comparer le resultat de filtrage avec le resultat précédant de FASTQC.

## Téléchargé des données à partir d'un lien
### Traiement des données en single end


1. importer les données NGS en utilsant ces liens
````
https://zenodo.org/record/3977236/files/female_oral2.fastq-4143.gz
````
<img width="1478" height="570" alt="image" src="https://github.com/user-attachments/assets/4a332378-d3b8-4855-a6db-c22137fc716e" />


2. Lancer `FastQC` et inspecter le resultats
3. Comment vous trouvez ces données en terme de qualité?
4. Utiliser blastn/VecScreen pour voir quels sont les dapatateurs trouvés
5. Utiliser cutadapt pour enlever les adaptateurs et améliorer les sequences brutes ( [adapter list](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314))
6. Quel pourcentage de lectures contient un adaptateur ?
7. Quel pourcentage de lectures a été tronqué en raison d'une mauvaise qualité ?
8. Quel pourcentage de lectures a été supprimé car trop courtes ?

Nous pouvons examiner nos données tronquées avec  FastQC et faire une comparaison des resultats.
### Traitement des données en paired end
1. Telecharger ces données
````
https://zenodo.org/record/61771/files/GSM461178_untreat_paired_subset_1.fastq
https://zenodo.org/record/61771/files/GSM461178_untreat_paired_subset_2.fastq
````
2.  Lancer FASTQC
3.  Que pensez vous des données?
4.  Lancer cutadapat ou fastp pour améliorer la qualité.
5.  Combien de pb enlevées a cause de la qualité basse?
6.  Combien de séquences courtes supprimées?
---

# 🧬Assemblage génomique
L'assemblage génomique consiste à reconstruire la séquence complète d'un génome à partir de millions de fragments courts d'ADN (les reads) obtenus par séquençage. C'est un puzzle informatique géant à résoudre sans modèle de départ
Dans cette étape nous allons faire plusier assemblage génomique sur différents jeux de données.



## Cas d'étude 1 :Assemblage d’un génome bactérien MRSA avec des données Illumina MiSeq

> Source : galaxy training
---

<a id="objectifs"></a>
## 1. Objectifs

À la fin de cet atelier, les participants pourront :

- importer des reads Illumina paired-end dans Galaxy ;
- organiser les données en collection paired-end ;
- évaluer la qualité des reads ;
- nettoyer les reads avec `fastp` ;
- assembler un génome bactérien avec `Shovill` ;
- évaluer l’assemblage avec `QUAST`.

---

<a id="contexte-biologique"></a>
## 2. Contexte biologique

Cet atelier utilise des données de séquençage Illumina MiSeq provenant d’une souche de **Staphylococcus aureus résistant à la méthicilline** (**MRSA**).

Le MRSA est un pathogène important impliqué dans les infections hospitalières.  
Assembler son génome permet ensuite de réaliser plusieurs analyses :

- recherche de gènes de résistance aux antibiotiques ;
- comparaison entre souches ;
- typage génomique ;
- étude épidémiologique.

---

<a id="preparation-historique"></a>
## 3. Préparation de l’historique Galaxy

Avant de commencer, créer un nouvel historique Galaxy.

### Étapes

1. Ouvrir Galaxy.
2. Créer un nouvel historique.
3. Renommer l’historique, par exemple :

```text
MRSA_Illumina_Assembly
```

<details>
<summary>Pourquoi créer un nouvel historique ?</summary>

Un historique séparé permet d’organiser l’analyse, d’éviter de mélanger les fichiers et de rendre le travail plus reproductible.

</details>

---

<a id="importation-donnees"></a>
## 4. Importation des données

Importer les deux fichiers FASTQ paired-end suivants :

```text
https://zenodo.org/record/10669812/files/DRR187559_1.fastqsanger.bz2
https://zenodo.org/record/10669812/files/DRR187559_2.fastqsanger.bz2
```

### Étapes dans Galaxy

1. Cliquer sur **Upload Data**.
2. Choisir **Paste/Fetch Data**.
3. Coller les deux liens.
4. Cliquer sur **Start**.
5. Fermer la fenêtre.

### Renommer les datasets

Renommer les fichiers en :

```text
DRR187559_1
DRR187559_2
```

### Créer une collection paired-end

1. Sélectionner les deux fichiers.
2. Choisir **Build List of Dataset Pairs**.
3. Vérifier l’appariement `_1` / `_2`.
4. Nommer la collection :

```text
Paired Reads
```

5. Ajouter le tag :

```text
#unfiltered
```

---

## Question — Que contient un fichier FASTQ ?

<details>
<summary>Afficher la réponse</summary>

Un fichier FASTQ contient quatre lignes par read :

1. un identifiant commençant par `@` ;
2. la séquence nucléotidique ;
3. une ligne séparatrice `+` ;
4. les scores de qualité encodés en ASCII.

Exemple :

```text
@read001
ATGCGTACGTAG
+
IIIIHHHHFFFF
```

</details>

---

<a id="controle-qualite"></a>
## 5. Contrôle qualité des reads

Le séquençage peut introduire plusieurs erreurs :

- mauvaise qualité en fin de read ;
- erreurs d’appel de bases ;
- adaptateurs résiduels ;
- reads trop courts ;
- biais de composition GC.

Avant l’assemblage, il est donc nécessaire de contrôler la qualité des données.

### Outil utilisé

Dans le tutoriel original, l’outil utilisé est :

```text
Falco
```

Falco est une alternative rapide à FastQC pour les données courtes de type Illumina.

### Paramètre Galaxy

Outil : **Falco**

```text
Raw read data from your current history: Paired Reads
```

### Résultats à examiner

Dans le rapport HTML, observer :

- Per base sequence quality ;
- Per sequence quality scores ;
- Per base sequence content ;
- GC content ;
- Adapter content.

---

## Interprétation rapide : Per base sequence quality

- Axe X : position de la base dans le read.
- Axe Y : score de qualité Phred.
- Vert : bonne qualité.
- Jaune : qualité moyenne.
- Rouge : mauvaise qualité.

<details>
<summary>Pourquoi la qualité baisse souvent vers la fin des reads Illumina ?</summary>

Au cours du séquençage Illumina, le signal fluorescent devient parfois moins net avec les cycles.  
Cela provoque souvent une baisse progressive de qualité vers la fin des reads.

</details>

---

<a id="nettoyage-fastp"></a>
## 6. Nettoyage avec fastp

Après le contrôle qualité, les reads peuvent être nettoyés avec `fastp`.

### Objectif

Le trimming permet de retirer :

- les bases de mauvaise qualité ;
- les adaptateurs ;
- les reads trop courts.

### Paramètres Galaxy

Outil : **fastp**

```text
Single-end or paired reads: Paired Collection
Select paired collection(s): Paired Reads
```

### Filtrage par longueur

```text
Length required: 30
```

### Coupe selon la qualité

```text
Cut by quality in front (5’): Yes
Cut by quality in tail (3’): Yes
Cutting window size: 4
Cutting mean quality: 20
```

### Rapport

```text
Output JSON report: Yes
```

Après fastp :

- supprimer le tag `#unfiltered` ;
- ajouter le tag `#filtered`.

---

## Question — Quel est l’effet attendu du trimming ?

<details>
<summary>Afficher la réponse</summary>

Le trimming peut :

- réduire légèrement la longueur moyenne des reads ;
- augmenter la proportion de bases Q20 et Q30 ;
- supprimer les adaptateurs ;
- améliorer la qualité globale des données.

Le contenu GC ne doit normalement pas changer fortement.

</details>

---

<a id="assemblage-shovill"></a>
## 7. Assemblage avec Shovill

Une fois les reads nettoyés, on peut assembler le génome.

### Pourquoi Shovill ?

`Shovill` est un assembleur basé sur SPAdes, optimisé pour les petits génomes bactériens.

### Paramètres Galaxy

Outil : **Shovill**

```text
Input reads type: Paired Collection
Paired collection: sortie paired-end de fastp
Estimated genome size: 2914567
```

### Sorties principales

Shovill produit :

1. un fichier log ;
2. un fichier FASTA contenant les contigs ;
3. un graphe d’assemblage.

---

## Question — Qu’est-ce qu’un contig ?

<details>
<summary>Afficher la réponse</summary>

Un contig est une séquence continue reconstruite à partir de reads qui se chevauchent.

Avec des reads courts Illumina, un génome bactérien est souvent assemblé en plusieurs contigs.

</details>

---

<a id="evaluation-quast"></a>
## 8. Évaluation avec QUAST

`QUAST` permet d’évaluer la qualité d’un assemblage génomique.

### Paramètres Galaxy

Outil : **QUAST**

```text
Assembly mode: Co-assembly
Contigs/scaffolds file: sortie Contigs de Shovill
```

### Métriques importantes

Observer dans le rapport QUAST :

- nombre de contigs ;
- longueur totale de l’assemblage ;
- N50 ;
- taille du plus grand contig ;
- contenu GC.

---

## Question — Comment interpréter un rapport QUAST ?

<details>
<summary>Afficher la réponse</summary>

Un bon assemblage bactérien devrait avoir :

- une longueur totale proche de la taille attendue du génome ;
- un nombre limité de contigs ;
- un N50 élevé ;
- un contenu GC cohérent avec l’espèce.

Cependant, avec des reads courts, il est fréquent de ne pas obtenir un génome fermé en un seul contig.

</details>

---

<a id="conclusion"></a>
## 9. Conclusion

Workflow réalisé :

```text
Importation FASTQ
        ↓
Contrôle qualité
        ↓
Nettoyage avec fastp
        ↓
Assemblage avec Shovill
        ↓
Évaluation avec QUAST
```

### Points clés

- Le contrôle qualité est indispensable avant l’assemblage.
- `fastp` améliore la qualité des reads.
- `Shovill` est adapté aux petits génomes bactériens.
- Les reads courts peuvent donner un assemblage fragmenté.
- `QUAST` permet d’évaluer quantitativement l’assemblage.

---


## Cas d'étude 2 :Assemblage de plusieurs génomes bacteriens



<a id="objectifs"></a>
## 1. Objectifs

À la fin de cet atelier, les participants pourront :

- importer des reads Illumina paired-end à partir d'un historique Galaxy ;
- évaluer la qualité des reads ;
- nettoyer les reads avec `fastp` ;
- assembler un génome bactérien avec `Shovill` ;
- évaluer l’assemblage avec `QUAST`.
- Annoter les génomes avec ``Prokka``
- Visualiser les genome sous ``jbrowse``

---

<a id="contexte-biologique"></a>
## 2. Contexte biologique

Cet atelier utilise des données de séquençage Illumina provenant d’une souche de **Salmonella galinarum** agent pathogen de la typhose aviaire.
 
Assembler son génome permet ensuite de réaliser plusieurs analyses :
- recherche de gènes de résistance aux antibiotiques ;
- comparaison entre souches ;
- typage génomique ;
- étude épidémiologique.

---

<a id="preparation-historique"></a>
## 3. Préparation de l’historique Galaxy

Avant de commencer, créer un nouvel historique Galaxy.


---

<a id="importation-donnees"></a>
## 4. Importation des données
### Étapes

1. Ouvrir Galaxy.
2. ouvirie le lien ci dissous de l'historique partagé
   ``https://usegalaxy.eu/u/solay/h/inh-workshop-2026``
2. Importer l'historique entier

<img width="905" height="500" alt="image" src="https://github.com/user-attachments/assets/3485978e-26ae-4bb5-b253-33c32a958cd2" />

3. Vous pouvez renomer l'historique


## 5. Contrôle qualité des reads

Le séquençage peut introduire plusieurs erreurs :

- mauvaise qualité en fin de read ;
- erreurs d’appel de bases ;
- adaptateurs résiduels ;
- reads trop courts ;
- biais de composition GC.

Avant l’assemblage, il est donc nécessaire de contrôler la qualité des données.

### Outil utilisé

Dans le tutoriel original, l’outil utilisé est :

```text
Falco ou FastQC
```

Falco est une alternative rapide à FastQC pour les données courtes de type Illumina.

### Paramètre Galaxy

Outil : **Falco**

```text
Raw read data from your current history: Paired Reads
```

### Résultats à examiner

Dans le rapport HTML, observer :

- Per base sequence quality ;
- Per sequence quality scores ;
- Per base sequence content ;
- GC content ;
- Adapter content.

---

## 6. Nettoyage avec fastp

Après le contrôle qualité, les reads peuvent être nettoyés avec `fastp`.

### Objectif

Le trimming permet de retirer :

- les bases de mauvaise qualité ;
- les adaptateurs ;
- les reads trop courts.

### Paramètres Galaxy

Outil : **fastp**

```text
Single-end or paired reads: Paired Collection
Select paired collection(s): Paired Reads
```
### Filtrage par longueur

```text
Length required: 30
```
### Coupe selon la qualité

### Rapport

```text
Output JSON report: Yes
```

## Question — Quel est l’effet attendu du trimming ?

<details>
<summary>Afficher la réponse</summary>

Le trimming peut :

- réduire légèrement la longueur moyenne des reads ;
- augmenter la proportion de bases Q20 et Q30 ;
- supprimer les adaptateurs ;
- améliorer la qualité globale des données.

Le contenu GC ne doit normalement pas changer fortement.

</details>

---

<a id="assemblage-shovill"></a>
## 7. Assemblage avec Shovill

Une fois les reads nettoyés, on peut assembler le génome.

### Pourquoi Shovill ?

`Shovill` est un assembleur basé sur SPAdes, optimisé pour les petits génomes bactériens.


### Sorties principales

Shovill produit :

1. un fichier log ;
2. un fichier FASTA contenant les contigs ;
3. un graphe d’assemblage.

---

<a id="evaluation-quast"></a>
## 8. Évaluation avec QUAST

`QUAST` permet d’évaluer la qualité d’un assemblage génomique.

### Métriques importantes

Observer dans le rapport QUAST :

- nombre de contigs ;
- longueur totale de l’assemblage ;
- N50 ;
- taille du plus grand contig ;
- contenu GC.

---

## Question — Comment interpréter un rapport QUAST ?

<details>
<summary>Afficher la réponse</summary>

Un bon assemblage bactérien devrait avoir :

- une longueur totale proche de la taille attendue du génome ;
- un nombre limité de contigs ;
- un N50 élevé ;
- un contenu GC cohérent avec l’espèce.

Cependant, avec des reads courts, il est fréquent de ne pas obtenir un génome fermé en un seul contig.

</details>

---
## 9. Annotation avec prokka
L'annotation génomique est le processus qui consiste à identifier les éléments fonctionnels d'une séquence d'ADN après son assemblage. Si l'assemblage permet d'obtenir les "pages du livre" du génome, l'annotation en écrit l'index et le dictionnaire pour en comprendre le sens.

Prokka est le logiciel standard de bioinformatique utilisé pour l'annotation rapide et automatisée de génomes de procaryotes (bactéries et archées). Développé par Torsten Seemann, son principal avantage est sa vitesse : il peut annoter un génome bactérien complet en moins de 10 minutes.

<img width="308" height="434" alt="image" src="https://github.com/user-attachments/assets/d4fd26ee-ab5a-4445-b862-880e424ee237" />

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



