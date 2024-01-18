from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
import sys

def fetch_sequence(accession):
    Entrez.email = "yx546@cornell.edu"
    try:
        handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=accession)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
    except Exception as e:
        sys.stderr.write(f"Error fetching data: {e}\n")
        return None

def compare_scores(cpp_file_name, python_score):
    try:
        with open(cpp_file_name, 'r') as file:
            for line in file:
                if "score" in line.lower():
                    cpp_score = float(line.split(':')[-1].strip())
                    if cpp_score == python_score:
                        print("success")
                    else:
                        print("fail")
    except FileNotFoundError:
        print(f"File {cpp_file_name} not found. Please ensure the C++ program has been run.")


def smith_waterman_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty):
    aligner = PairwiseAligner()
    aligner.mode = 'local'  # local alignment
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_penalty
    aligner.extend_gap_score = gap_penalty

    alignments = aligner.align(seq1, seq2)
    return alignments[0]

accession1 = "GQ457487"
seq1 = fetch_sequence(accession1)
accession2 = "GQ117044"
seq2 = fetch_sequence(accession2)

if seq1 is not None and seq2 is not None:

    # Save the sequences to files
    with open("sequence1.txt", "w") as file1, open("sequence2.txt", "w") as file2:
        file1.write(seq1)
        file2.write(seq2)

    # Run the Smith-Waterman algorithm and print the result
    alignment = smith_waterman_alignment(seq1, seq2, 3, -3, -2)
    print(alignment.score)

    # Compare the scores of the Python and C++ programs
    compare_scores("cpp_output.txt", alignment.score)
    
else:
    print("Error fetching sequences. Please check the error message above.")
