def save_fasta(sequences, filename, prefix="seq"):
    with open(filename, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">{prefix}{i}\n")
            f.write(seq + "\n")


def read_fasta(filename):
    headers = []
    sequences = []

    with open(filename, "r") as f:
        current_header = None
        current_sequence = []

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    headers.append(current_header)
                    sequences.append("".join(current_sequence))

                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line.upper())

        if current_header is not None:
            headers.append(current_header)
            sequences.append("".join(current_sequence))

    return headers, sequences

def save_fasta_with_headers(headers, sequences, filename):
    with open(filename, "w") as f:
        for header, seq in zip(headers, sequences):
            f.write(f">{header}\n")
            f.write(seq + "\n")


def make_meme_safe_headers(headers, prefix="seq"):
    safe_headers = []
    for i, _ in enumerate(headers, start=1):
        safe_headers.append(f"{prefix}_{i}")
    return safe_headers