Pour organisation et reflexion : https://annuel2.framapad.org/p/CMamIs*2

# PROJET PYTHON CMI
## I - Consignes
### Enoncé
2 rendu, un oral et un code commenté.
10 min d'oral le projet.
Faire un planning


L’algorithme de programmation dynamique a été utilisé dès 1970 pour aligner les séquences biologiques. Il existe plusieurs variantes : alignement global (Needleman & Wunsh, 1970), alignement semi-global (Sellers, 1974) et alignement local (Waterman, 1981).

L’objectif de ce projet est de programmer cet algorithme fondateur dans sa version la plus simple c’est-à-dire celle où la fonction de pénalisation des indels est une constante. Un exemple de calcul est donné dans la feuille excel PDmatrices.xls.

Récursion, basée sur le principe d'optimalité : 
hTps://fr.wikipedia.org/wiki/Algorithme_de_Needleman-Wunsch

Vous trouverez également sur ceTe page wikipédia le pseudo-code du calcul de la matrice.
```
for i=0 to length(A)-1
    F(i, 0) ← d*i
for j=0 to length(B)-1
    F(0,j) ← d*j 
for i=1 to length(A)-1
	for j = 1 to length(B)-1 
	{
	 Choice1 ← F(i-1,j-1) + S(A(i), B(j))
	 Choice2 ← F(i-1, j) + d
	 Choice3 ← F(i, j-1) + d
	 F(i, j) ← max(Choice1, Choice2, Choice3)
	}
```
Vous rendrez le code qui permet de calcul la matrice des coûts et celle des chemins (cf feuille excel) ainsi que celui qui permet à parrir de la matrice des chemins de reconstruire l’alignement.
- Etape 1 : écrire le programme pour aligner des séquences nucléiques
- Etape 2 : améliorer la sortie (cf annexe) et l’enregistrer dans un fichier de sortie
- Etape 3 : écrire un programme pour aligner des séquences protéiques
- Etape 4 : écrire un programme qui détecte automatiquement le type de séquences et envoie le bon programme avec les bons paramètres.

PS : les paramètres du programme peuvent être fournis dans un fichier de configuration 
```
#======================================= # # Aligned_sequences: 2
# 1: INS1_MOUSE
# 2: INS2_MOUSE
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 110
# Identity: 97/110 (88.2%)
# Similarity: 99/110 (90.0%)
# Gaps: 2/110 ( 1.8%)
# Score: 501.5
#
#
#=======================================
INS1_MOUSE   1 MALLVHFLPLLALLALWEPKPTQAFVKQHLCGPHLVEALYLVCGERGFFY 50
               |||.:.||||||||.|||..||||||||||||.|||||||||||||||||
INS2_MOUSE   1 MALWMRFLPLLALLFLWESHPTQAFVKQHLCGSHLVEALYLVCGERGFFY 50
INS1_MOUSE  51 TPKSRREVEDPQVEQLELGGSP--GDLQTLALEVARQKRGIVDQCCTSIC 98
               ||.||||||||||.||||||.|  |||||||||||:||||||||||||||
INS2_MOUSE  51 TPMSRREVEDPQVAQLELGGGPGAGDLQTLALEVAQQKRGIVDQCCTSIC 100 INS1_MOUSE  99 SLYQLENYCN 108
               ||||||||||
INS2_MOUSE 101 SLYQLENYCN 110
#---------------------------------------
#---------------------------------------
```
ANNEXE: extrait de la sortie produite par le logiciel NEEDLE de la suite logicielle EMBOSS.

## II - Notre plan de programmation

### About our code

#### Auteurs
Anais Bosc-Bierne, Blandine Hottekiet-Genetier, Enora Cieslak, Helena Drude, Mathis Bourgoin, Paul Lemonier

Emails: aboscbierne@etud.univ-angers.fr , bhott@etud.univ-angers.fr , enciesl@etud.univ-angers.fr ,helena.drude@etud.univ-angers.fr , mbourgoin@etud.univ-angers.fr , paul.lemonnier@etud.univ-angers.fr

#### Licences

### MODULE (Fichier fonctions)

#### Fonction 1
```
Calcul matrice coûts
+ Affichage dans la console #optionnel (booléen "debug"= FALSE par défaut)
```

#### Fonction 2
```
Calcul alignement
+ Affichage dans la console
+ Sortie ouptut dans Fichier sortie #optionnel (booléen "expert" = FALSE par défaut)
```

### Programme principal (fichier ouvert par l'utilisateur)

1- Ouverture fichiers fasta
    -> seq1
    -> seq2

2- Comptage des nucléotides
    => pour vérifier que la taille de la séquence ne dépasse pas 20 nt

3- Question à l'utilisateur pour qu'il choisisse entre les 3 modes (débuggage =1, courant=2, expert=3) :
    ```
    -> si 1) :
    appel fonction 1 avec debug=TRUE
    mat.cost<-(x) #avec x = le return de la fonction 1
    appel fonction 2
    
    -> si 2) :
    appel fonction 1 avec debug=FALSE
    mat.cost<-(x) #avec x = le return de la fonction 1
    appel fonction 2

    -> si 3) :
    appel fonction 1 avec debug=FALSE
    mat.cost<-(x) #avec x = le return de la fonction 1
    appel fonction 2 avec expert=TRUE
    ```
    



