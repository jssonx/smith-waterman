from Bio.Align import PairwiseAligner

def compare_scores(cpp_file_name, python_score):
    try:
        with open(cpp_file_name, 'r') as file:
            for line in file:
                if "score" in line.lower():
                    cpp_score = float(line.split(':')[-1].strip())
                    if cpp_score == python_score:
                        print("Success: The score match.")
                    else:
                        print("Fail: The score does not match.")
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

def read_sequence_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read().strip()
    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return None

seq1_path = "./data/GQ457487.txt"
seq2_path = "./data/GQ117044.txt"
seq1 = read_sequence_from_file(seq1_path)
seq2 = read_sequence_from_file(seq2_path)

if seq1 is not None and seq2 is not None:

    alignment = smith_waterman_alignment(seq1, seq2, 3, -3, -2)
    print(f"Alignment score: {alignment.score}")

    compare_scores("cpp_output.txt", alignment.score)
else:
    print("Error reading sequences. Please check the files in the data directory.")
