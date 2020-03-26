# GeneLink test code.

import os

# Assume that the project file structure is:
# dna/
#    data/     # There can be multiple data directories 
#    pyprog/

dataDir = "data/"

# List the .fasta files in a specified directory: ../dataDir
# It may be useful to include other extensions that are used
# for fasta files (e.g., fns)
def ld(verbose=True, dataDir = dataDir):
   files = os.listdir("../" + dataDir)
   fileList = []
   for file in files:
      if file.endswith(".fasta"):  
         fileList.append(file)
         if verbose:
            print (file)
   return filesList

# Aminio acid code letter to name
def a2n(a = 'F'):
   if a == 'F':
      return "phenylalanine"
   elif a == 'L':
      return "leucine"
   elif a == 'I':
      return "isoleucine"
   elif a == 'M':
      return "methionine"
   elif a == 'V':
      return "valine"
   elif a == 'S':
      return "serine"
   elif a == 'P':
      return "proline"
   elif a == 'T':
      return "threonine"
   elif a == 'A':
      return "alanine"
   elif a == 'Y':
      return "tyrosine"
   elif a == 'H':
      return "histidine"
   elif a == 'Q':
      return "glutamine"
   elif a == 'N':
      return "asparagine"
   elif a == 'K':
      return "lysine"
   elif a == 'D':
      return "aspartic acid"
   elif a == 'E':
      return "glutamic acid"
   elif a == 'C':
      return "cysteine"
   elif a == 'W':
      return "tryptophan"
   elif a == 'R':
      return "arginine"
   elif a == 'S':
      return "serine"
   elif a == 'G':
      return "glycine"

# Transform a codon to an amino acid code letter.  Return '' for stop codons.
# Note that 'ATG' returns 'M' (Methionine) rather than 'Start' since it is
# assumed that the function will be called only with codons occuring in introns.
# According to [Science  323(5911), 259-61. 2009 Jan 9] and possibly other references,
# the stop codon TGA can code for an amino acid in certain cases.  The code
# here does not account for such insances.
def c2a(c = 'TTT'):
    if len(c) < 3:
       return ''
    if c in ['TTT', 'TTC']:
        return 'F'
    elif c in ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']:
        return 'L'
    elif c in ['ATT', 'ATC', 'ATA']:
        return 'I'
    elif c in ['ATG']: # Note that ATG is also the start codon
        return 'M'
    elif c in ['GTT', 'GTC', 'GTA', 'GTG']:
        return 'V'
    elif c in ['TCT', 'TCC', 'TCA', 'TCG']:
        return 'S'
    elif c in ['CCT', 'CCC', 'CCA', 'CCG']:
        return 'P'
    elif c in ['ACT', 'ACC', 'ACA', 'ACG']:
        return 'T'
    elif c in ['GCT', 'GCC', 'GCA', 'GCG']:
        return 'A'
    elif c in ['TAT', 'TAC']:
        return 'Y'
    elif c in ['CAT', 'CAC']:
        return 'H'
    elif c in ['CAA', 'CAG']:
        return 'Q'
    elif c in ['AAT', 'AAC']:
        return 'N'
    elif c in ['AAA', 'AAG']:
        return 'K'
    elif c in['GAT', 'GAC']:
        return 'D'
    elif c in ['GAA', 'GAG']:
        return 'E'
    elif c in ['TGT', 'TGC']:
        return 'C'
    elif c in ['TGG']:
        return 'W'
    elif c in ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']:
        return 'R'
    elif c in ['AGT', 'AGC']:
        return 'S'
    elif c in ['GGT', 'GGC', 'GGA', 'GGG']:
        return 'G'
    elif c in ['TAA', 'TAG']:    # Ochre stop
        return 'stop'
    elif c in ['TAG']:           # Amber Stop
        return 'stop'
    elif c in ['TGA']:           # Opal Stop (see above comments about this codon)
        return 'stop'
    return 'error'

# Transform a DNA sequence into a sequence of amino acids.  The input
# should be a string of ACGT (or acgt) but may contain other characters
# (which are ignored).
def seq2c2(sequence, verbose = False):
    seq    = sequence['sequence']
    n      = len(seq)    
    name   = sequence['name']
    header = sequence['header']

    amino_acids   = []
    triple        = ''

    # Stop and start indices are stored as the index of first molecule in the codon.
    # For example, each would be stored as index 7 in the following:
    # ATG    (start codon)
    # TGA    (a stop codon)
    # 789
    starts        = []        # list of indices the start codons ATG
    stops         = []        # list of indices the stop codons
    
    caseSensitive = ''        # This will hold case sensitive encoding (cap = coding)
    introns       = []        # This will contain a list of introns (each is a string)
    n_introns     = 0
    n_expons      = 0
    error_char_positions = [] # List of positions that are not A, C, G or T
    error_chars   = []        # List of the erroneous characters
    error_codes   = []
    
    # State of the calculation
    WAITING_FOR_START   = 0
    READING_AMINO_ACIDS = 1
    state               = WAITING_FOR_START
    start_index         = 0
    i                   = 0  # Current position in the sequence.
    
    while i < len(seq):
       if state == WAITING_FOR_START:
          # Prepare to start a new amino acid
          amino_acid = ''
          
          if (i + 3) > len(seq):   # Can't read a full codon.
             # We are in a noncoding region.  Finish off the sequence in lower case.
             caseSensitive = caseSensitive + seq[i:len(seq)].lower()
             i = len(seq)
          else:                     # Can read a full codon.
             triple = seq[i:(i+3)]  # Read next codon.
             if triple == 'ATG':
                #start_index   = i
                starts.append(i)
                caseSensitive = caseSensitive + 'atg'  #  start seq is noncoding
                i = i + 3   # We must have read three bases since the triple was ATG
                state = READING_AMINO_ACIDS
             else: # Codon was not a start sequence

                # Add the character to the sequence even if it isn't among ACGTacgt.
                # Note that some files already use upper/lower case.
                caseSensitive = caseSensitive + seq[i].lower()

                # Some files have characters (errors?) not among ACGTacgt.
                if seq[i] not in 'ACGTacgt': 
                   error_char_positions.append(i)
                   error_chars.append(seq[i])
                   
                i = i + 1   # Read noncoding bases one at a time
                
       elif state == READING_AMINO_ACIDS:
          if (i + 3) >= len(seq): # In coding region but can't read a full codon.
             # We are in a coding region.  Finish off the sequence in upper case.
             caseSensitive = caseSensitive + seq[i:n].upper()
             i = len(seq)

             # We can't read another amino acid since there aren't three bases left.  
             # Add what we have to the list of amino acids but also add an error code.
             amino_acids.append(amino_acid)
             error_codes.append('Incomplete final codon of length ' + str(len(seq) - i) +
                                " in coding region;")
             
          else:                      # Can read a full codon.
             codon    = seq[i:(i+3)]
             i        = i + 3             
             aa       = c2a(codon)

             if aa == 'stop':        # We read a stop.
                stops.append(i-3)
                caseSensitive = caseSensitive + codon.lower()
                amino_acids.append(amino_acid)
                state = WAITING_FOR_START
             else:                   # Read either an amino acid or an error.
                if aa in 'FLIMVSPTAYHQNKDECWRSG':  # We read an amino acid.
                   amino_acid    = amino_acid + aa
                   
                elif aa == 'error':  # Codon was not an amino acid.
                   # Add a blank to the amino acid if we read three bases in a coding 
                   # region but at least one base is not among ACGTacgt.
                   amino_acid    = amino_acid + ' '
                   error_codes.append("Read codon " + codon + " at position " + str(i-1) +
                                      " in coding region " + str(len(introns)+1) + "; ")

                # Keep upper case (coding) even if one or more bases was not among ACGTacgt.
                caseSensitive = caseSensitive + codon.upper()
       
    # Compute lengths of coding and noncoding sections
    coding_length    = 0
    noncoding_length = 0
    for x in caseSensitive:
       if x == 'A' or x == 'C' or x == 'G' or x == 'T':
          coding_length = coding_length + 1
       elif x == 'a' or x == 'c' or x == 'g' or x == 't':
          noncoding_length = noncoding_length + 1
       else:
          # If other characters occur in the sequence (such as N), assume they
          # are noncoding.
          noncoding_length = noncoding_length + 1 
    
    coding_fraction    = coding_length / len(seq)
    noncoding_fraction = noncoding_length / len(seq)

    # compute number of introns
    introns = len(starts)
    
    # compute number of exons
    exons = introns + 1

    return {"name":name,
            "header":sequence['header'],
            "filename":sequence['filename'],
            "amino_acids":amino_acids,
            "starts":starts,
            "stops":stops,
            "sequence":caseSensitive,
            "introns":introns,
            "introns":introns,
            "coding_length":coding_length,
            "coding_fraction":coding_fraction,
            "exons":exons,
            "noncoding_length":noncoding_length,
            "noncoding_fraction":noncoding_fraction,
            "error_chars":error_chars,
            "error_char_positions":error_char_positions,
            "error_codes":error_codes}


# Process all data files.
def pf(dataDir = dataDir):
   files = ld(dataDir = dataDir)
   for file in files():
      s = readseq2(file, dataDir = dataDir)
      #convert to 

def myflatten(list_of_strings = []):
   u = ''
   for x in list_of_strings:
      u = u + x
   return u

# Compute statistics about the sequence.
# Inputs:
#   a processed sequence
# Outputs:
#   a dictionary with the following keys
#      length: length of the sequence
#      i: coding percent of the sequence
#      e: noncoding percent of the sequence
#      A: percent of coding region that is A
#      C: percent of coding region that is C
#      G: percent of coding region that is G
#      T: percent of coding region that is T
#      a: percent of noncoding region that is a
#      c: percent of noncoding region that is c
#      g: percent of noncoding region that is g
#      t: percent of noncoding region that is t
#      C1: percent of coding region that is amino acid C
#      D1: percent of coding region that is amino acid D
#      E1: percent of coding region that is amino acid E
#      F1: percent of coding region that is amino acid F
#      G1: percent of coding region that is amino acid G
#      H1: percent of coding region that is amino acid H
#      I1: percent of coding region that is amino acid I
#      K1: percent of coding region that is amino acid K
#      L1: percent of coding region that is amino acid L
#      M1: percent of coding region that is amino acid M
#      N1: percent of coding region that is amino acid N
#      P1: percent of coding region that is amino acid P
#      Q1: percent of coding region that is amino acid Q
#      R1: percent of coding region that is amino acid R
#      S1: percent of coding region that is amino acid S
#      T1: percent of coding region that is amino acid T
#      V1: percent of coding region that is amino acid V
#      Y1: percent of coding region that is amino acid Y
#      W1: percent of coding region that is amino acid W
#
def stat(seq, verbose=False):
   results = {}
   if verbose:
      print(seq['filename'])
   results['length'] =  len(seq['sequence'])
   results['i']      = seq['coding_fraction']
   results['e']      = seq['noncoding_fraction']
   results['A']      = 0.0
   results['C']      = 0.0
   results['G']      = 0.0
   results['T']      = 0.0
   results['a']      = 0.0
   results['c']      = 0.0
   results['g']      = 0.0
   results['t']      = 0.0
   results['A1']     = 0.0
   results['C1']     = 0.0
   results['D1']     = 0.0
   results['E1']     = 0.0
   results['F1']     = 0.0
   results['G1']     = 0.0
   results['H1']     = 0.0
   results['I1']     = 0.0
   results['K1']     = 0.0
   results['L1']     = 0.0
   results['M1']     = 0.0
   results['N1']     = 0.0
   results['P1']     = 0.0
   results['Q1']     = 0.0
   results['R1']     = 0.0
   results['S1']     = 0.0
   results['T1']     = 0.0
   results['V1']     = 0.0
   results['Y1']     = 0.0
   results['W1']     = 0.0
   for x in seq['sequence']:
       if (not (x == 'n') and not (x in 'rwyslkm')):
         results[x] = results[x] + 1
   results['A'] = results['A'] / seq['coding_length']
   results['C'] = results['C'] / seq['coding_length']
   results['G'] = results['G'] / seq['coding_length']
   results['T'] = results['T'] / seq['coding_length']
   results['a'] = results['a'] / seq['noncoding_length']
   results['c'] = results['c'] / seq['noncoding_length']
   results['g'] = results['g'] / seq['noncoding_length']
   results['t'] = results['t'] / seq['noncoding_length']
   for x in myflatten(seq['amino_acids']):
       key = x + "1"
       results[key] = results[key] + (1 / len(seq['amino_acids']))
   return (results)

def test3(filename, dataDir = dataDir):
   seq  = readseq2(filename, dataDir = dataDir)
   seq1 = seq2c2(seq)
   coding_length    = 0
   noncoding_length = 0   
   for x in seq1['cs']:
      if x == 'A' or x == 'C' or x == 'G' or x == 'T':
         coding_length = coding_length + 1
      if x == 'a' or x == 'c' or x == 'g' or x == 't':
         noncoding_length = noncoding_length + 1
   print ([coding_length, seq1['coding_length'], coding_length - seq1['coding_length']])
   print ([noncoding_length, seq1['noncoding_length'], noncoding_length - seq1['noncoding_length']])

def test2(filename, dataDir = dataDir):
   seq  = readseq2(filename, dataDir = dataDir)
   seq1 = seq2c2(seq)
   st   = stat(seq1)
   print([st['A'] + st['C'] + st['G'] + st['T']])   
   print([st['a'] + st['c'] + st['g'] + st['t']])   
    
def test1(filename, dataDir = dataDir):
   seq  = readseq2(filename, dataDir = dataDir)
   seq1 = seq2c2(seq)
   A  = C  = G  = T  = 0
   A1 = C1 = G1 = T1 = 0
   for x in seq['sequence']:
      if x == 'A':
        A = A + 1
      elif x == 'C':
        C = C + 1
      elif x == 'G':
        G = G + 1
      elif x == 'T':
        T = T + 1
   for x in seq1['cs']:
      if x == 'A' or x == 'a':
         A1 = A1 + 1
      if x == 'C' or x == 'c':
         C1 = C1 + 1
      if x == 'G' or x == 'g':
         G1 = G1 + 1
      if x == 'T' or x == 't':
         T1 = T1 + 1
   print ([A, A1, A - A1])
   print ([C, C1, C - C1])
   print ([G, G1, G - G1])
   print ([T, T1, T - T1]) 
   print ([len(seq['sequence']), len(seq1['cs']), len(seq['sequence']) - len(seq1['cs'])])
   print ([len(seq['sequence']), seq1['coding_length'] + seq1['noncoding_length'], len(seq['sequence']) - seq1['coding_length'] - seq1['noncoding_length']])

# Extract a short name and full name from the header of a fasta file.
def getName(header):
   n1        = header.find(' ')  # Find the first space character in the first line
   n         = len(header)
   shortname = header[1:n1]      # Chop off the initial > character.
   name      = header[(n1+1):n]
   return shortname, name
   
# Read A DNA sequence from a FASTA formatted file.  The first
# line read should be of the form
# >string0 string1
# It is also assumed that the data file is in dataDir where
# dataDir is defined at the top of this file.
# Inputs:
#      filename: the name of the data file
#      verbose : if True, print the name of the file that was read
#      
# Outputs:
#      a dictionary with two keys:
#         filename: the name of the file that was read
#         header  : the complete first line of the file
#         name    : the name computed from getName(header)
#         sequence: the genome sequence
def readseq2(filename, verbose=False, dataDir = dataDir):
   path = "../" + dataDir + "/" + filename
   seq       = ""
   f         = open(path, 'r')
   header    = f.readline().strip() # Take the entire first line as header

   shortname, name = getName(header)
   
   for line in f:
      seq = seq + line.rstrip() # Remove white space from the end of the line
   f.close()
   
   if verbose:
      print(filename + " read")
      
   return {"sequence":seq, "name":name, "filename":filename, "header":header}


# Read all the fasta files in a specified file list.
# If no file list is specified, all files in dataDir
# are read.  
# Inputs:
#    fileList: a list of file names to be found in dataDir
# Outputs:
#    a list of the sequences read from the files
def readAll(fileList = [], dataDir = dataDir):
   os.chdir("../" + dataDir)
   if len(fileList) == 0:
      fileList = ld(verbose=False, dataDir=dataDir)
   else:
      pass # Read only files from the list
   sequences = []
   for file in fileList:
      sequences.append(readseq2(file, dataDir=dataDir))
   os.chdir("../pyprog")
   return sequences

import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from sklearn.decomposition import PCA
def runpca(fileList = [], dataDir=dataDir):
   features = ['i', 'e', 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'C1', 'D1', 'E1', 'F1', 'G1',
               'H1', 'I1', 'K1', 'L1', 'M1', 'N1', 'P1', 'Q1', 'R1', 'S1', 'T1', 'V1', 'W1']
   sequences = readAll(fileList, dataDir=dataDir)
   sequences2 = []
   for s in sequences:
      print (s['name'])
      sequences2.append(stat(seq2c2(s)))
   x = [[0]* len(features)  for i in range(len(sequences2))]
   for i in range(len(sequences2)):
      for j in range(len(features)):
         x[i][j] = sequences2[i][features[j]]
   X = np.array(x)
   pca = PCA(n_components = 2)
   principalComponents = pca.fit_transform(X)
   principalDF = pd.DataFrame(data = principalComponents, columns = ['pc1', 'pc2'])
   principalDF.to_csv("principalDF.csv")
   #pca.fit(X)
   print(type(principalDF))
   plt.scatter(principalDF['pc1'], principalDF['pc2'])
   plt.show()   
   #return pca

#from sklearn.cluster import KMeans
# Apply k-means clustering baded on a specified feature list.
# [incomplete]
def km(features = ['i'], fileList = [], dataDir = dataDir):
   sequences = readAll(fileList, dataDir=dataDir)
   sequences2 = []
   for s in sequences:
      sequences2.append(stat(seq2c2(s)))
   x = [[0]* len(features)  for i in range(len(sequences2))]
   for i in range(len(sequences2)):
      for j in range(len(features)):
         x[i][j] = sequences2[i][features[j]]
   return x

# Plot two dimensions of features.  That is, read all the sequence 
# data in the data directory, compute all the statistics then plot
# two specified dimensions.
def p2d(dims = [0, 1], features = ['i', 'C1'], fileList = [], dataDir = dataDir):
   if len(dims) < 2 or len(features) < 2:
      print ("I need two features to make a plot.  Goodbye.")
      return
   if len(dims) > 2:
      print ("I can only make a plot with two features.  Going with the first two.")
   dims2 = dims[0:2]
   points = km(features, fileList, dataDir=dataDir)
   x = [points[i][dims2[0]] for i in range(len(points))]
   y = [points[i][dims2[1]] for i in range(len(points))]
   plt.scatter(x, y)
   plt.xlabel(features[dims2[0]])
   plt.ylabel(features[dims2[1]])
   plt.show()
      
