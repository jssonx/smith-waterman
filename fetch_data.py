from Bio import Entrez, SeqIO
import sys
import os

def fetch_sequence(accession):
    Entrez.email = "yx546@cornell.edu"
    try:
        handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=accession)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequence = str(record.seq)

        if not os.path.exists("./data"):
            os.makedirs("./data")

        with open(f"./data/{accession}.txt", "w") as file:
            file.write(sequence)

        print(f"Sequence {accession} saved to ./data/{accession}.txt")
        return sequence
    except Exception as e:
        sys.stderr.write(f"Error fetching data: {e}\n")
        return None

fetch_sequence("GQ457487")
fetch_sequence("GQ117044")