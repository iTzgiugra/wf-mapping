#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Input files
params.reference = "reference.fasta"
params.sample = "sample.fastq"

// Define the workflow
workflow {

    // Process to align reads
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

    // Process to convert SAM to BAM
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

    // Process to generate statistics
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

    // Process to calculate alignment percentage
    process calculatePercentage {
        label 'python'
        container 'my-python-image' // Ensure the Docker image is correctly built and available

        input:
        path "alignment_stats.txt"

        output:
        path "alignment_percentage.txt"

        script:
        """
        cp /usr/src/app/calculate_percentage.py .
        python calculate_percentage.py alignment_stats.txt alignment_percentage.txt
        """
    }

    // Workflow steps
    alignedSam = alignReads(params.reference, params.sample)
    alignedBam = convertToBam(alignedSam)
    alignmentStats = generateStats(alignedBam)
    calculatePercentage(alignmentStats)
}
