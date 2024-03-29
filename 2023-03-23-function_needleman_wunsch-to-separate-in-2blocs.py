## ---------------------------
##
## Script name: projet
##
## Purpose of script: Projet CMI Programmation_Landes
##
## Version: 1.0
##
## Author: Anais Bosc-Bierne, Blandine Hottekiet-Genetier, Enora Cieslak, Helena Drude, Mathis Bourgoin, Paul Lemonier
## Email: aboscbierne@etud.univ-angers.fr , bhott@etud.univ-angers.fr , enciesl@etud.univ-angers.fr ,
## helena.drude@etud.univ-angers.fr , mbourgoin@etud.univ-angers.fr , paul.lemonnier@etud.univ-angers.fr
##
## ---------------------------

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

def needleman_wunsch(seq1, seq2, iseq1="", iseq2="", match_score=4, mismatch_penalty=-1, gap_penalty=-2):
    ### 1) FOR THE FUNCTION "CALCULATION OF THE COSTS MATRIX"  ##
    
    # Initialisation of the similarities matrix == costs matrix #
    #############################################################
    
    ## Calculating the lengths of the two sequences
    n = len(seq1)
    m = len(seq2)
    
    ## Initialisation of the matrix variable
      #as a 2D list of zeros with n+1 rows and m+1 columns
      #+1 to add the gaps' scores in the matrix
    matrix = [[0 for j in range(m+1)] for i in range(n+1)]
    
    ## Filling of the 1st column : with gaps' scores
      #In other words : the 1st sequence is aligned to gaps in the 2nd seq
    for i in range(1, n+1):
        #Increasing gap penalties as the length of the gap increases
        matrix[i][0] = i * gap_penalty
    ## Filling of the 1st line : with gaps' scores
      #In other words : the 2nd sequence is aligned to gaps in the 1st seq
    for j in range(1, m+1):
        matrix[0][j] = j * gap_penalty

    # Filling of the costs matrix #
    ###############################
    
    ## Reading of the 1st sequence, base by base (lines)
    for i in range(1, n+1):
      ## And Reading of the 2nd sequence in parallel, base by base (columns)
        for j in range(1, m+1):
            ## Calculation of the 3 possible scores (depending on the way of alignment)
            
            #Score if diagonal way == identity (match) or substitution (mismatch)
              #previous code : diagonal_score = matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
              #here I would rather put seq1[i] than seq1[i-1] and seq2[j] rather than seq2[j-1], as corrected below
            diagonal_score = matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            
            #Score if horiztontal way == indel (gap)
              #with deletion//gap in the seq1, so we move by column but the line remains the same
            left_score = matrix[i][j-1] + gap_penalty
            
            #Score if vertical way == indel (gap)
              #with deletion//gap in the seq2, so we move by line but the column remains the same
            up_score = matrix[i-1][j] + gap_penalty
            
            ## Registration of the highest score (sum of costs) in the cell of the costs matrix 
            matrix[i][j] = max(diagonal_score, left_score, up_score)
            
    ### 2) TO ADD : FUNCTION "TO PRINT THE COSTS MATRIX""  ##
    
    ### 3) FOR THE FUNCTION "CALCULATION OF THE ALIGNMENT""  #####
    
    # Alignment : Trace back through the costs matrix to find it #
    #############################################################
   #'''This part of the code is responsible for tracing back through the matrix
   #to find the optimal alignment between the two input sequences.
   #The algorithm starts at the bottom right corner of the matrix
   #(corresponding to the end of both sequences) and works its way
   #towards the top left corner (corresponding to the beginning of both sequences).'''
   
    ## Initialisation of the new sequences, for the alignment
    #as empty strings
    align1 = ""
    align2 = ""
    #Starting i and j to the ends of their sequences, to begin the traceback
    i = n
    j = m
    
    #this will wok until none of the sequence is finished
    while i > 0 and j > 0:
        ## Comparison of the 3 possible scores to find the good way of alignment
            
        #Score if diagonal way == identity (match) or substitution (mismatch)
          #previous code : matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
          #here I would rather put seq1[i] than seq1[i-1] and seq2[j] rather than seq2[j-1], as corrected below
        if matrix[i][j] == matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            #adding the element of each sequence in its alignement version
              #thanks to the order of the concatenation : seq[] + align and not the reverse (align + seq[]),
              #each iteration will fill the alignment version of the sequence
              #from the end of the sequence to the beginning, so in the same order as we read it
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            #we move backward in each sequence ; given that i -=1 <=> i=i-1
            i -= 1
            j -= 1
            
        #Score if horizontal way == indel (gap)
          #with deletion//gap in the seq1 (because the line i remains i)...
        elif matrix[i][j] == matrix[i][j-1] + gap_penalty:
            #adding the gap in the alignment version of seq1
            align1 = "-" + align1
            #adding the element of seq2 in its alignment version
            align2 = seq2[j-1] + align2
            #...so we move (backward) in the column but the line remains the same
            j -= 1
            
        #Last option : Score if vertical way == indel (gap)
          #with deletion//gap in the seq2...
        else:
             #adding the element of seq1 in its alignment version
            align1 = seq1[i-1] + align1
            #adding the gap in the alignment version of seq2
            align2 = "-" + align2
            #...so we move (backward) in the line but the column remains the same
            i -= 1
            
    #this will wok if seq1 is finished (i=0) but seq2 is not (j>0)
      #until seq2 is also finished
    while j > 0:
        #filling of seq1 with gaps
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j -= 1

    #this will wok if seq2 is finished (j=0) but seq1 is not (i>0)
      #until seq1 is also finished
    while i > 0:
        align1 = seq1[i-1] + align1
        #filling of seq2 with gaps
        align2 = "-" + align2
        i -= 1
        
    '''To sum up the result of this part of the code : it's a traceback algorithm
    to determine the optimal alignment % the two input sequences
    based on the costs matrix that was previously filled.
    The results are 2 strings which represent the aligned sequences,
    where dashes indicate gaps in the alignment.'''

    # Calculation of the identities, similarities, and gaps #
    #########################################################
    
    ## Length of the alignment
      # total length = length of align1 = length of align2
    length = len(align1)
    
    ## Counting the identities
      #iterating over the aligned sequences : "zip(,)"
      ##counting a amount (sum) of positions where the elements (=characters) are the same
    identity = sum(1 for a, b in zip(align1, align2) if a == b)
    
    ## Counting the similarities
      #counting a amount (sum) of the number of positions where the characters are the same
      #or at least one is a gap
    similarity = sum(1 for a, b in zip(align1, align2) if a == b or a == '-' or b == '-')
    
    ## Counting the gaps
    gaps = sum(1 for a, b in zip(align1, align2) if a == '-' or b == '-')

    # Calculation of the score
    ##########################
    "At the end of the loop, the variable 'score' holds the final score of the alignment."
    
    ## Initialisation of the score
    score = 0
    
    ## Reading the alignment sequences from their beginning to their end
    for i in range(length):
        #Score if gap
        if align1[i] == '-' or align2[i] == '-':
            #given that score += x <=> score = score + x
            score += gap_penalty
        #Score if identity (match)
        elif align1[i] == align2[i]:
            score += match_score
        #Score if substitution (mismatch)
        else:
            score += mismatch_penalty
    
    ### 4) FOR THE FUNCTION "FINAL PRINTING" ##############
    ### THE ALIGNMENT NEEDS TO END IN AN EXTERNAL FILE* ###
    ### *IN ADDITION OF THE OUTPUT PRINTED IN THE CONSOLE #
    # Output : formatting like EMBOSS Needle ##############
    #######################################################
    
    output = ""
    output += f"\nPairwise alignment of :\n"
    
    ## Recovering the variables defining the sequences ; ex: "iseq1" : what is it ??
      #this variable does not seem to be defined in this code
      #would this be the information on the first line of the fasta file of the sequence ?
    output += f"Sequence 1 : {iseq1}\n"
    output += f"Sequence 2 : {iseq2}\n"
    
    output += f"Using Needleman-Wunsch algorithm\n"
    
    ## Recovering the parameters (costs) of the alignment
      #from the variables of the code
    output += f"Match score of {match_score}\n"
    output += f"Mismatch penalty of {mismatch_penalty}\n"
    output += f"Gap penalty of {gap_penalty}\n"
    
    ## Recovering the characteristics of the obtained alignment
    output += f"\nLength: {length}\n"
      #Calculating the percentage of identities, similarities and gaps
        #the function parameter .2f limits the number of digits after the comma
        #(for a decimal number) - I think
    output += f"Identity: {identity}/{length} ({100 * identity / length:.2f}%)\n"
    output += f"Similarity: {similarity}/{length} ({100 * similarity / length:.2f}%)\n"
    output += f"Gaps: {gaps}/{length} ({100 * gaps / length:.2f}%)\n"
    output += f"Score: {score}\n\n"
    
    ## Initialsation of the nucleotides counters
    #for the seq1 (align1)
    num1_nt = 0
    #for the seq2 (align2)
    num2_nt = 0
    
    ## Printing of the aligned sequences
    #with no more than 50 characters* by row, so i is going from 50 to 50
    for i in range(0, len(align1), 50):
        # * here we see below align1[i:i+50] to represent this limit of 50 nt
        #{num_nt1} is for printing the variable which counts the number of nt
        output += f"Sequence 1 : {num1_nt}\t{align1[i:i+50]}\t"
        
        #to count the number of nt, we have to substact the gaps...
          # - which we count thanks to ".count()" -
          #...from the total number of characters
        num1_nt = num1_nt + len(align1[i:i+50]) - align1[i:i+50].count("-")
        
        #total number of nt, for seq1
        output += f"{num1_nt}\n"
        
        #to put spaces below "Sequence 1 :" to allow the alignment between the sequences and the "|" line
        output += f" "*13
        
        #tabulation to allow the alignment between the sequences
        output += f"\t"
        
        #printing | to indicate the identities
        output += f"".join("|" if a == b else " " for a, b in zip(align1[i:i+50], align2[i:i+50])) + "\n"
        
        #printing the alignment version of seq2 (=align2)
        output += f"Sequence 2 : {num2_nt}\t{align2[i:i+50]}\t"
        
        num2_nt = num2_nt + len(align2[i:i+50]) - align2[i:i+50].count("-")
        
        #total number of nt, for seq2
        output += f"{num2_nt}\n"        
    output += f"\n"
    return output

seq1 = "ATTCAAGCTGA"
seq2 = "AACTTGCGTGA"
fin = needleman_wunsch(seq1, seq2, iseq1="", iseq2="", match_score=4, mismatch_penalty=-1, gap_penalty=-2)
print(fin)
