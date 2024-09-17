# Alignment Workflow

## Beschreibung

Dieser Workflow verwendet BWA zum Ausrichten von Reads an ein Referenzgenom, konvertiert die SAM-Dateien mit Samtools in BAM-Dateien, erstellt Alignment-Statistiken und visualisiert die Ergebnisse mit einem Python-Skript.

## Voraussetzungen

- Nextflow
- Docker
- BWA
- Samtools
- Python 3.x
- matplotlib, pandas, seaborn (Python-Pakete)

## Verwendung

1. Klone das Repository:
    ```bash
    git clone https://github.com/iTzgiugra/wf-mapping.git
    ```

2. Navigiere in das Verzeichnis:
    ```bash
    cd wf-mapping/wf-mapping
    ```

3. Führe den Workflow aus:
    ```bash
    nextflow run main.nf --reference /home/graziano/reference/reference.fasta --sample /home/graziano/probe/sample.fastq
    ```

## Lizenz

Dieses Projekt ist unter der [MIT-Lizenz](LICENSE) lizenziert.
