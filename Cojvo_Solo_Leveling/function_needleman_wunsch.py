def needleman_wunsch(seq1, seq2, iseq1="", iseq2="", match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    # Initialize the matrix
    n = len(seq1)
    m = len(seq2)
    matrix = [[0 for j in range(m+1)] for i in range(n+1)]
    for i in range(1, n+1):
        matrix[i][0] = i * gap_penalty
    for j in range(1, m+1):
        matrix[0][j] = j * gap_penalty

    # Fill in the matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            diagonal_score = matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            left_score = matrix[i][j-1] + gap_penalty
            up_score = matrix[i-1][j] + gap_penalty
            matrix[i][j] = max(diagonal_score, left_score, up_score)

    # Trace back through the matrix to find the alignment
    align1 = ""
    align2 = ""
    i = n
    j = m
    while i > 0 and j > 0:
        if matrix[i][j] == matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i][j-1] + gap_penalty:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
        else:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1

    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j -= 1

    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = "-" + align2
        i -= 1

    # Calculate the identity, similarity, and gaps
    length = len(align1)
    identity = sum(1 for a, b in zip(align1, align2) if a == b)
    similarity = sum(1 for a, b in zip(align1, align2) if a == b or a == '-' or b == '-')
    gaps = sum(1 for a, b in zip(align1, align2) if a == '-' or b == '-')

    # Calculate the score
    score = 0
    for i in range(length):
        if align1[i] == '-' or align2[i] == '-':
            score += gap_penalty
        elif align1[i] == align2[i]:
            score += match_score
        else:
            score += mismatch_penalty
    
    # Format the output like EMBOSS Needle
    output = ""
    output += f"\nPairwise alignment of :\n"
    output += f"Sequence 1 : {iseq1}\n"
    output += f"Sequence 2 : {iseq2}\n"
    output += f"Using Needleman-Wunsch algorithm\n"
    output += f"Match score of {match_score}\n"
    output += f"Mismatch penalty of {mismatch_penalty}\n"
    output += f"Gap penalty of {gap_penalty}\n"
    output += f"\nLength: {length}\n"
    output += f"Identity: {identity}/{length} ({100 * identity / length:.2f}%)\n"
    output += f"Similarity: {similarity}/{length} ({100 * similarity / length:.2f}%)\n"
    output += f"Gaps: {gaps}/{length} ({100 * gaps / length:.2f}%)\n"
    output += f"Score: {score}\n\n"
    
    num_gaps1 = 0
    num_gaps2 = 0
    for i in range(0, len(align1), 50):
        output += f"Sequence 1 : {num_gaps1}\t{align1[i:i+50]}\t"
        num_gaps1 = num_gaps1 + align1[i:i+50] - align1[i:i+50].count("-")
        output += f"{num_gaps1}\n"
        output += f" "*13
        output += f"\t"
        output += f"".join("|" if a == b else " " for a, b in zip(align1[i:i+50], align2[i:i+50])) + "\n"
        output += f"Sequence 2 : {num_gaps2}\t{align2[i:i+50]}\t"
        num_gaps2 = num_gaps2 + align2[i:i+50] - align2[i:i+50].count("-")
        output += f"{num_gaps2}\n"        
    output += f"\n"
    return output
