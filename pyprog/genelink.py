
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
      return "aspartic"
   elif a == 'E':
      return "glutamic Acid"
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

def c2a(c = 'TTT'):
    if c in ['TTT', 'TTC']:
        return 'F'
    elif c in ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']:
        return 'L'
    elif c in ['ATT', 'ATC', 'ATA', 'ATG']:
        return 'I'
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
    elif c in ['TAA', 'TAG']:    # Stop
        return ''
    elif c in ['TGA']:           # Stop
        return ''
      
def seq2c():
    amino_acids = []
    # read 3 characters at a time.
    # while c not = 'ATG', continue.
    # amino_acids.append(c2a(c))
    return amino_acids
