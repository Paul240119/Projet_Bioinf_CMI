import tkinter as tk
from tkinter import filedialog

def choose_file():
    # create a Tkinter root window, but hide it
    root = tk.Tk()
    root.withdraw()

    # prompt the user to select a file and get the file path
    file_path = filedialog.askopenfilename()

    # read the contents of the selected file
    with open(file_path, 'r') as file:
        file = file.read()

    # return the text content of the file
    return file

# define function to parse the 1st input file
def open_seq1():
    
    ## Initialisation of the variables...
    #...which will be used by other functions (function 2)...
    #...so we create them as global
    global seq1
    global iseq1
    
    # call choose_file_and_read_text() to get the text content of the selected file
    file = choose_file()

    # split the text into lines and get the first line
    lines = file.split('\n')
    first_line = lines[0]

    # check if the first line starts with ">"
    if first_line.startswith('>'):
        print("Information : The first line starts with '>'")
        iseq1 = first_line
        seq1 = '\n'.join(lines[1:])
        print("iseq1:",iseq1)
        print("seq1:",seq1)
    else:
        print("Information : The first line does not start with '>'")
        iseq1 = ''
        seq1 = file
        print("iseq1:",iseq1)
        print("seq1:",seq1)

def open_seq2():
    
    ## Initialisation of the variables...
    #...which will be used by other functions (function 2)...
    #...so we create them as global
    global seq2
    global iseq2
    
    # call choose_file to get the text content of the selected file
    file = choose_file()

    # split the text into lines and get the first line
    lines = file.split('\n')
    first_line = lines[0]

    # check if the first line starts with ">"
    if first_line.startswith('>'):
        print("Information : The first line starts with '>'")
        iseq2 = first_line
        seq2 = '\n'.join(lines[1:])
        print("iseq2:",iseq2)
        print("seq2:",seq2)
    else:
        print("Information : The first line does not start with '>'")
        iseq2 = ''
        seq2 = file
        print("iseq2:",iseq2)
        print("seq2:",seq2)

def verif_seq(seq):
    allowed_chars = set(["A", "T", "C", "G","\n"])
    seq_length = len(seq)

    # check if all characters in seq are in allowed_chars set
    if set(seq) <= allowed_chars:
        print("Valid DNA sequence")

        # check if seq length is less than or equal to 20
        if seq_length <= 20:
            print("Length of seq is less than or equal to 20")
            propose_debug = True
            
        else:
            print("Length of seq is greater than 20")
            propose_debug = False

    else:
        print("Warning : Give a valid DNA sequence")
        propose_debug = False
        quit()
        
    return propose_debug

def mode_choice():
    
    #PREPARING THE WINDOW FOR INTERACTIVE QUESTION TO THE USER
    root = tk.Tk()
    root.withdraw()

    # create a new window
    mode_choice_window = tk.Toplevel(root)

    #set the window title
    mode_choice_window.title("Mode Choice")

    # set the question for the prompt
    prompt_label = tk.Label(mode_choice_window, text="Select a mode:")
    prompt_label.pack()

    #variable to store the mode which will be selected
    mode_var = tk.StringVar(mode_choice_window)

    #PREPARING THE RESULT OF THE FUNCTION
    #which is the mode choice of the user
    #previous code : de select_mode(choice)
    def select_mode():
        mode_var.set() #define choice depending on the button pressed (which is =mode_var)
        mode_choice_window.destroy()
    
    #PREPARING THE BUTTONS which may be displayed after
    debug_button = tk.Button(mode_choice_window, text="Debug Mode", 
                             command=lambda: select_mode('debug'))
    normal_button = tk.Button(mode_choice_window, text="Normal Mode",
                             command=lambda: select_mode('normal'))
    expert_button = tk.Button(mode_choice_window, text="Expert Mode",
                             command=lambda: select_mode('expert'))
    
    #CREATING THE BUTTONS
    #if the sequences are not too long (propose_debug=True for both)
    if verif1 and verif2:
        #displaying the buttons for debug mode, normal mode, and expert mode
        debug_button.pack()
        normal_button.pack()
        expert_button.pack()
    
    #if one of the sequence is too long
    else:
        #not displaying the button for debug mode
        normal_button.pack()
        expert_button.pack()
    
    #DISPLAYING THE WINDOW, WITH THE 2 OR 3 BUTTONS
    mode_choice_window.wait_window(mode_choice_window)  # wait for user to make selection
    
    #RETURNING the choice of the user ('debug', 'normal' or 'expert')
    return mode_var.get()

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
    #global matrix
    
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
              #/!\ there is a shift between i(matrix) and i(seq),
                #because of the "-" (indel) line and column of the matrix
                #it means that : element[i] in the matrix <=> element[i-1] in the sequence
            diagonal_score = matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            
            #Score if horiztontal way == indel (gap)
              #with deletion//gap in the seq1, so we move by column but the line remains the same
            left_score = matrix[i][j-1] + gap_penalty
            
            #Score if vertical way == indel (gap)
              #with deletion//gap in the seq2, so we move by line but the column remains the same
            up_score = matrix[i-1][j] + gap_penalty
            
            ## Registration of the highest score (sum of costs) in the cell of the costs matrix 
            matrix[i][j] = max(diagonal_score, left_score, up_score)

    # Optionnal* printing of the costs matrix #
    ##########################################

    #*If the user chooses the debug mode
    if debug==True :
        #Calculate the width of the sequence 2
        width=((len(seq2)+1)*6)+5
        
        #printing seq2
        print("{:^5s}{:^6s}".format("","-"),end="")
        for char in seq2:
            print("{:^6s}".format(char),end="")
        print("\n","*"*width,sep="")
        
        #initialisation of the counter
        i=0
        #Reading the matrix (which is a sort of list of lists), list by list
        #it means that here each "line" in our false dataframe is a list for python,
        #and the whole of these lists make our false dataframe (matrix)
        for list in matrix:
            
            # TITLE OF THE LINE #
            #For all the lists except the first one
            if i!=0:
                #printing the nucleotide n°i(matrix)-1 of seq1
                print("{:^3s}{:^3s}".format(seq1[i-1],"|"),end="")
                #preparing i in order to read the next list at next turn
                i=i+1
                
            #For the 1st list
            if i==0:
                #printing the - for indel (deletion in seq1)
                print("{:^3s}{:^3s}".format("-","|"),end="")
                #preparing i in order to read the next list at next turn
                i=i+1

            # FILLING OF THE LINE #
            #printing all the costs of this list...
            for element in list:
                #...which correspond to the possible alignments...
                #...of the element of seq1[i-1] with all the elements of seq2[]
                print("{:^3d}{:^3s}".format(element,"|"),end="")
            print("\n",end="")
            print("-"*width)

    return matrix

def alignment(matrix, seq1, seq2, expert =False,  iseq1="", iseq2="", match_score=4, mismatch_penalty=-1, gap_penalty=-2):
    
    # Alignment : Trace back through the costs matrix to find it #
    #############################################################
   #'''This part of the code is responsible for tracing back through the matrix
   #to find the optimal alignment between the two input sequences.
   #The algorithm starts at the bottom right corner of the matrix
   #(corresponding to the end of both sequences) and works its way
   #towards the top left corner (corresponding to the beginning of both sequences).'''
   
    ## Initialisation of the new sequences, for the alignment
    #as empty strings
    n = len(seq1)
    m = len(seq2)
    align1 = ""
    align2 = ""
    #Starting i and j to the ends of their sequences, to begin the traceback
    i = n
    j = m
    
    #this will work until none of the sequence is finished
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
            
    #this will work if seq1 is finished (i=0) but seq2 is not (j>0)
      #until seq2 is also finished
    while j > 0:
        #filling of seq1 with gaps
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j -= 1

    #this will work if seq2 is finished (j=0) but seq1 is not (i>0)
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
    
    ### FOR THE "FINAL PRINTING" ##############
    ### THE ALIGNMENT CAN END IN AN EXTERNAL FILE* ###
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
    
    if expert == False:
          print(output) #Printing the results in the console
    else:
       with open("Alignement.txt", "w") as filout:
           filout.write(output) #Creating a file with the alignement 


def Needleman() :
    
    #Import the sequences
    open_seq1()
    open_seq2()
    
    #verification of the sequences
    global verif1 #to use it in other functions (mode_choice)
    global verif2
    verif1=verif_seq(seq1)
    verif2=verif_seq(seq2)

    #Mode Choice of the user
    choice = mode_choice() 

    if choice == "debug" :
        #will return matrix and print it in the console
        matrix_cost(seq1, seq2, match_score=4, mismatch_penalty=-1, gap_penalty=-2, debug=True)
        #will make the alignment and print it in the console too
        alignment(matrix, seq1, seq2, expert =False,  iseq1="", iseq2="", match_score=4, mismatch_penalty=-1, gap_penalty=-2)
    
    elif choice == "expert" :
        #will return matrix (without printing it in the console)
        matrix_cost(seq1, seq2, match_score=4, mismatch_penalty=-1, gap_penalty=-2, debug=False)
        #will make the alignment and write it in an external file (without printing it in the console)
        alignment(matrix, seq1, seq2, expert=True,  iseq1="", iseq2="", match_score=4, mismatch_penalty=-1, gap_penalty=-2)
    
    #if choice == "normal"
    else :
        #will return matrix (without printing it in the console)
        matrix_cost(seq1, seq2, match_score=4, mismatch_penalty=-1, gap_penalty=-2, debug=False)
        #will make the alignment and print it in the console)
        alignment(matrix, seq1, seq2, expert=False,  iseq1="", iseq2="", match_score=4, mismatch_penalty=-1, gap_penalty=-2)

