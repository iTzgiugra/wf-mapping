import sys

def parse_flagstat(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if 'mapped (' in line:
                percentage = line.split('(')[1].split('%')[0].strip()
                return float(percentage)
    raise ValueError("Mapped percentage not found in the file.")

def main():
    if len(sys.argv) != 3:
        print("Usage: python calculate_percentage.py <stats_file> <output_file>")
        sys.exit(1)

    stats_file = sys.argv[1]
    output_file = sys.argv[2]

    percentage = parse_flagstat(stats_file)

    with open(output_file, 'w') as f:
        f.write(f"Alignment percentage: {percentage}%\n")

if __name__ == "__main__":
    main()
