import os

def extract(filename, dataDir):
    os.chdir("../" + dataDir)
    f = open(filename, 'r')
    outfile = open("tempfile.fasta", 'w')
    for line in f:
        if line[0] == ">":
            outfile.close()
            name = line.split(' ')[0]
            name = name[1:len(name)] + ".fasta"
            print (name)
            outfile = open(name, 'w')
        outfile.write(line)
    f.close()
    os.remove("tempfile.fasta")
    os.chdir("../pyprog")
    
