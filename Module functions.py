def matrix_cost(seq1, seq2, match_score=4, mismatch_penalty=-1, gap_penalty=-2, debug=False):
    # PARAMETERS OF MATRIX_COST #
        #seq1:sequence FASTA n°1
        #seq2:sequence FASTA n°2

        #match_score,mismatch_penalty,gap_penalty : parameters of the costs

        #debug : boolean to choose to print the matrix of costs or not (at the end) ;
            #default value (False)=does not print it
    
    # Initialisation of the similarities matrix == costs matrix #
    #############################################################

    ## Initialisation of the variables...
    #...which will be used by other functions (function 2)...
    #...so we create them as global
    global n
    global m
    global matrix
    
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

    # Optionnal printing of the costs matrix #
    ###############################

    if debug==True :
        width=((len(seq2)+1)*6)+5
        
        #printing seq2
        print("{:^5s}{:^6s}".format("","-"),end="")
        for char in seq2:
            print("{:^6s}".format(char),end="")
        print("\n","*"*width,sep="")
        
        i=0
        for list in matrix:
            if i!=0:
                #printing the nucleotide n°i of seq1
                print("{:^3s}{:^3s}".format(seq1[i-1],"|"),end="")
                i=i+1
            if i==0:
                #printing the - for indel
                print("{:^3s}{:^3s}".format("-","|"),end="")
                i=i+1
            for element in list:
                #printing the cost corresponding to the alignment
                print("{:^3d}{:^3s}".format(element,"|"),end="")
            print("\n",end="")
            print("-"*width)

    return matrix

## Test appel fonction ci-dessous
seq1="ATTCAAGCTGA"
seq2="AACTTGCGTGA"

matrix_cost(seq1, seq2, match_score=4, mismatch_penalty=-1, gap_penalty=-2,debug=True)








    

