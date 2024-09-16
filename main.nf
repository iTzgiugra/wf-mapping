#!/usr/bin/env nextflow

// Input files
params.reference = "reference.fasta"
params.sample = "sample.fastq"

// Define the workflow
workflow {

    process alignReads {
        label 'bwa'
        input:
        path reference
        path sample

        output:
        path "aligned.sam"

        script:
        """
        bwa index $reference
        bwa mem $reference $sample > aligned.sam
        """
    }

    process convertToBam {
        label 'samtools'
        input:
        path "aligned.sam"

        output:
        path "aligned.bam"

        script:
        """
        samtools view -S -b aligned.sam > aligned.bam
        """
    }

    process generateStats {
        label 'samtools'
        input:
        path "aligned.bam"

        output:
        path "alignment_stats.txt"

        script:
        """
        samtools flagstat aligned.bam > alignment_stats.txt
        """
    }

    process plotResults {
        label 'python'
        input:
        path "alignment_stats.txt"

        output:
        path "alignment_report.png"

        script:
        """
        python3 bin/plot_alignment.py alignment_stats.txt alignment_report.png
        """
    }

    // Workflow steps
    alignReads(params.reference, params.sample)
    convertToBam()
    generateStats()
    plotResults()
}

