#!/usr/bin/env python

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match_score: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError

# a simple function to read the name and sequence from a file
# The file is expected to have just one contig/sequence. This function
# checks the assumption and complains if it is not the case.
def read_single_contig_fasta(filename):
    names = []
    sequences = []
    with open(filename, 'r') as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T"]:
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]

def smith_waterman(seq1, seq2, match, mismatch, gapopen, gapextend):
  # Initialize the matrix for dynamic programming
  matrix = []
  # Iterate over the length of seq2 to create rows
  for _ in range(len(seq2)+1):
      # Create a new row with length equal to the length of seq1, initialized to 0
      row = [0] * (len(seq1)+1)
      # Append the row to the matrix
      matrix.append(row)

  # Initialize the matrix for dynamic programming
  gap_sizes = []
  # Iterate over the length of seq2 to create rows
  for _ in range(len(seq2)+1):
      # Create a new row with length equal to the length of seq1, initialized to 0
      row = [0] * (len(seq1)+1)
      # Append the row to the matrix
      gap_sizes.append(row)

  # Initialize the matrix for dynamic programming
  traceback = []
  # Iterate over the length of seq2 to create rows
  for _ in range(len(seq2)+1):
      # Create a new row with length equal to the length of seq1, initialized to 0
      row = [0] * (len(seq1)+1)
      # Append the row to the matrix
      traceback.append(row)


  # Initialize variables to keep track of the maximum score and its position
  max_score = 0
  max_i = 0
  max_j = 0

  # Fill in the matrix using dynamic programming
  for i in range(1, len(seq2) + 1):
      for j in range(1, len(seq1) + 1):
          match_score = matrix[i - 1][j - 1] + (match if seq2[i - 1] == seq1[j - 1] else - mismatch)
          delete = (matrix[i - 1][j] - gapopen)
          insert = (matrix[i][j - 1] - gapopen)
          matrix[i][j] = max(0, match_score, delete, insert)
  
          # Update traceback matrix
          if matrix[i][j] == match_score:
              traceback[i][j] = 1  # Diagonal
          elif matrix[i][j] == insert:
              traceback[i][j] = 2  # left
          elif matrix[i][j] == delete:
              traceback[i][j] = 3  # up 

          # Update gap size matrix
          if traceback[i][j] == 2:  # left
              gap_sizes[i][j] = gap_sizes[1][j-1] + 1
          elif traceback[i][j] == 3:  # up
              gap_sizes[i][j] = gap_sizes[i-1][j] + 1
          else:
              gap_sizes[i][j] = 0

          # Update max score and its position          
          if matrix[i][j] > max_score:
              max_score = matrix[i][j]
              max_i = i
              max_j = j  
  
  print("scoring Matrix:")
  for row in matrix:
    print(row)
  
  
  # Print gap sizes matrix
  print("Gap Sizes Matrix:")
  for row in gap_sizes:
      print(row)

  # Calculate and print gap penalty for each position
  print("\nGap Penalty Calculation:")
  for i in range(1, len(seq2) + 1):
      for j in range(1, len(seq1) + 1):
          gap_penalty = gapopen + gap_sizes[i][j] * gapextend
          print(f"Position ({i}, {j}): Gap Penalty = {gap_penalty}")
  
  
  
  
  
  
  
  # Traceback to find the aligned sequences
  alnseq1 = ''
  alnseq2 = ''

   
  while traceback[max_i][max_j] != 0:
    if traceback[max_i][max_j] == 1:  # Diagonal
        alnseq1 = seq1[max_j - 1] + alnseq1
        alnseq2 = seq2[max_i - 1] + alnseq2
        max_i -= 1
        max_j -= 1
    elif traceback[max_i][max_j] == 2:  # Left
        alnseq1 = seq1[max_j - 1] + alnseq1
        alnseq2 = '-' + alnseq2
        max_j -= 1
    elif traceback[max_i][max_j] == 3:  # Up
        alnseq1 = '-' + alnseq1
        alnseq2 = seq2[max_i - 1] + alnseq2
        max_i -= 1  
   
  return max_score, alnseq1, alnseq2
    

def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and 
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(seq1, seq2, 
                                  match, mismatch, gapopen, gapextend)
    
    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:],
                     "hm:x:g:e:",
                     ["help", "match=", "mismatch=", "gapopen=", "gapextend="])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)
