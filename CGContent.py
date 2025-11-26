
import Toolkit


#data = """"""

res = {}
def getContent(data):
    for seq in data.split(">"):
        if seq != "":
            _, value = seq.split("_")
            number = value[:4]
            dna = value[4:]
            dna = dna.replace("\n", "")
            counts = Toolkit.count_nucleotides(dna)
            C = counts.get("C")
            G = counts.get("G")
            CG = C+G
            content = (CG / len(dna)) * 100
            res[number] = content
    num = max(res, key=res.get)
    print("Rosalind_" + num)
    print(res.get(num))





with open("rosalind_gc.txt", "r") as f:
    data = f.read()
getContent(data)