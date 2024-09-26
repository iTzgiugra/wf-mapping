from Bio import SeqIO
import sys

def find_matching_percentage(probe, reference):
    max_match = 0
    probe_len = len(probe)

    # Iteriere durch die Referenz und suche nach der Probe
    for i in range(len(reference) - probe_len + 1):
        match_count = 0
        # Prüfe jede Base der Probe auf Übereinstimmung
        for j in range(probe_len):
            if reference[i + j] == probe[j]:
                match_count += 1
            else:
                break  # Beende, wenn die Sequenz nicht mehr übereinstimmt
        max_match = max(max_match, match_count)

    # Berechne den Prozentsatz der Übereinstimmung
    matching_percentage = (max_match / probe_len) * 100
    return matching_percentage

def read_fastq(fastq_file):
    # Lese die Sequenz aus der FASTQ-Datei
    for record in SeqIO.parse(fastq_file, "fastq"):
        return str(record.seq)  # Nimmt die erste Sequenz aus der Datei

def read_fasta(fasta_file):
    # Lese die Sequenz aus der FASTA-Datei
    for record in SeqIO.parse(fasta_file, "fasta"):
        return str(record.seq)  # Nimmt die erste Sequenz aus der Datei

def main():
    if len(sys.argv) != 3:
        print("Usage: python calculate_percentage.py <probe_fastq> <reference_fasta>")
        sys.exit(1)

    probe_file = sys.argv[1]
    reference_file = sys.argv[2]

    # Lese die Sequenzen aus den Dateien
    probe_seq = read_fastq(probe_file)
    reference_seq = read_fasta(reference_file)

    # Berechne den Prozentsatz der Übereinstimmung
    percentage = find_matching_percentage(probe_seq, reference_seq)
    
    print(f"Matching percentage: {percentage:.2f}%")


if __name__ == "__main__":
    main()
