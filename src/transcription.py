#pour lancer ce petit script, vous devez passer la commande suivant dans le terminal
# python3 src/transcription.py data/sequence.fasta

import sys

sequence = ""

sequence_file = sys.argv[1]

with open(sequence_file,"r") as fh:
        for line in fh:
            if not line.startswith(">"):
                sequence = str(line)


cadre1 = []
for j in range(0,len(sequence),3):
    codon = sequence[j:j+3]
    cadre1.append(codon)


cadre2 = []
cadre2.append(sequence[0:1])
for j in range(1,len(sequence),3):
    codon = sequence[j:j+3]
    cadre2.append(codon)

cadre3 = []
cadre3.append(sequence[0:2])
for j in range(2,len(sequence),3):
    codon = sequence[j:j+3]
    cadre3.append(codon)


d = {
'TTT': 'F', 'TCT' : 'S' , 'TAT' : 'Y' , 'TGT' : 'C',
'TTC': 'F', 'TCC' : 'S' , 'TAC' : 'Y' , 'TGC' : 'C',
'TTA': 'L', 'TCA' : 'S' , 'TAA' : '*' , 'TGA' : '*',
'TTG': 'L', 'TCG' : 'S' , 'TAG' : '*' , 'TGG' : 'W',
'CTT': 'L', 'CCT' : 'P' , 'CAT' : 'H' , 'CGT' : 'R',
'CTC': 'L', 'CCC' : 'P' , 'CAC' : 'H' , 'CGC' : 'R',
'CTA': 'L', 'CCA' : 'P' , 'CAA' : 'Q' , 'CGA' : 'R',
'CTG': 'L', 'CCG' : 'P' , 'CAG' : 'Q' , 'CGG' : 'R',
'ATT': 'I', 'ACT' : 'T' , 'AAT' : 'N' , 'AGT' : 'S',
'ATC': 'I', 'ACC' : 'T' , 'AAC' : 'N' , 'AGC' : 'S',
'ATA': 'I', 'ACA' : 'T' , 'AAA' : 'K' , 'AGA' : 'R',
'ATG': 'M', 'ACG' : 'T' , 'AAG' : 'K' , 'AGG' : 'R',
'GTT': 'V', 'GCT' : 'A' , 'GAT' : 'D' , 'GGT' : 'G',
'GTC': 'V', 'GCC' : 'A' , 'GAC' : 'D' , 'GGC' : 'G',
'GTA': 'V', 'GCA' : 'A' , 'GAA' : 'E' , 'GGA' : 'G',
'GTG': 'V', 'GCG' : 'A' , 'GAG' : 'E' , 'GGG' : 'G',
}

def translate(cadre:list):
    word = ""
    for i in range(len(cadre)):
        if len(cadre[i]) > 2 :
            AA = d[cadre[i]]
            word = word + AA        
    return word

traduction1 = str(translate(cadre1))
traduction2 = str(translate(cadre2))
traduction3 = str(translate(cadre3))

print(traduction1)
print(traduction2)
print(traduction3)

def find_start_codon(traduction, which):
    for i in range(len(traduction)-1):
        word = ""
        if traduction[i] == "M":
            j = i+1
            while j < len(traduction):
                if traduction[j] == "*":
                    word = traduction[i:j]
                    break
                else:
                    j+=1 
        if len(word) > 30:
            print(word)
            return word
            
            
    

find_start_codon(traduction1,1)
find_start_codon(traduction2,2)
find_start_codon(traduction3,3)


