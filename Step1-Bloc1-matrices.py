#Code pour programmer le Bloc1 // Partie 1 de l'Etape 1
#=======================================
#
# Aligned_sequences: 2
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

#ouvre le fichier de la sequence 1 et celui de la sequence 2
#ouvre/cree un fichier matrice dans lequel ecrire
with open("C:/Users/bland/Downloads/Cours_230122/CMI_S2/Approfondissement programmation et modelisation/Projet/fichier_sequence1.txt","r") as fichier_sequence1,\
open("C:/Users/bland/Downloads/Cours_230122/CMI_S2/Approfondissement programmation et modelisation/Projet/fichier_sequence2.txt","r") as fichier_sequence2,\
open("C:/Users/bland/Downloads/Cours_230122/CMI_S2/Approfondissement programmation et modelisation/Projet/fichier_matrice.txt","w") as fichier_matrice :
#Pour chaque sequence lue dans le fichier "fichier_sequence1",
#on affiche son identifiant et son nom puis les 10 premiers nucleotides de sa sequence :
    sequence1 = ""
    for ligne in fichier_sequence1:
        if ligne.startswith(">"):
            sequence1 = ligne[1:].split()[0]
            sequence1_dict[sequence1] = ""
#utilise strip() pour supprimer le caractere de saut de ligne.
        else:
            sequence1[sequence1] += line.strip()
    for id in sequence1_dict:
        print(id)
        print(sequence1_dict[:10])
