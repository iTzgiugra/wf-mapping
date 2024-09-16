#!/usr/bin/env nextflow

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
        path alignedSam

        output:
        path "aligned.bam"

        script:
        """
        samtools view -S -b $alignedSam > aligned.bam
        """
    }

    // Process to generate statistics
    process generateStats {
        label 'samtools'
        input:
        path alignedBam

        output:
        path "alignment_stats.txt"

        script:
        """
        samtools flagstat $alignedBam > alignment_stats.txt
        """
    }

    // Process to plot results
    process plotResults {
        label 'python'
        input:
        path alignmentStats

        output:
        path "alignment_report.png"

        script:
        """
        python3 bin/plot_alignment.py $alignmentStats alignment_report.png
        """
    }

    // Workflow steps
    alignedSam = alignReads(params.reference, params.sample)
    alignedBam = convertToBam(alignedSam)
    alignmentStats = generateStats(alignedBam)
    plotResults(alignmentStats)
}
