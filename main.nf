#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    qcflow
    =================================
    Usage:
    abaaab
    Mandatory Arguments:
      --fastq              description
      --references         description
      --krakendb
    Options:
    Outputs:
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.fastq = false
params.references = false
params.krakendb = false
params.forward_adapter = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
params.reverse_adapter = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATC"

if ( params.fastq ) {
    fastqPairs = Channel.fromFilePairs(
        params.fastq,
        checkIfExists: true,
        size: 2,
        type: "file",
        flat: true
    )
} else {
    log.info "Hey I need some fastq files to look at please."
    exit 1
}


if ( params.references ) {
    references = Channel.fromPath(
        params.references,
        checkIfExists: true,
        type: "file"
    ) 
}

if ( !params.krakendb ) {

    process downloadKrakenDB {
        label "download"
        label "kraken"
        label "blast"

        storeDir "${params.outdir}/databases"

        output:
        file "krakendb" into krakenDB

        """
        mkdir -p krakendb

        kraken2-build --download-taxonomy --db krakendb

        kraken2-build --threads ${task.cpus} --download-library bacteria --db krakendb
        kraken2-build --threads ${task.cpus} --download-library archaea --db krakendb
        kraken2-build --threads ${task.cpus} --download-library viral --db krakendb
        kraken2-build --threads ${task.cpus} --download-library UniVec_Core --db krakendb
        kraken2-build --threads ${task.cpus} --download-library fungi --db krakendb

        kraken2-build \
            --threads ${task.cpus} \
            --build \
            --db krakendb \
            --kmer-len 35 \
            --minimizer-len 31 \
            --minimizer-spaces 6 \

        kraken2-build --clean
        """
    }

} else if ( file(params.krakendb).exists() ) {
    krakenDB = Channel.fromPath( params.krakendb )
} else {
    exit 1, "You specified a Kraken database, but it doesn't exist."
}


fastqPairs.into {
    fastqPairs4PretrimFastQC;
    fastqPairs4Cutadapt;
    fastqPairs4Bowtie;
    fastqPairs4Kraken;
}


process pretrimFastQC {
    label "java"
    label "fastqc"
    tag { fq.simpleName }

    publishDir "${params.outdir}/read_stats"

    input:
    set val(bName), val(direction), file(fq) from fastqPairs4PretrimFastQC
        .flatMap {bName, fwd, rev -> [[bName, "forward", fwd], [bName, "reverse", rev]]}

    output:
    set val(bName), val(direction), file("${fq.simpleName}_pretrim_fastqc.html"),
        file("${fq.simpleName}_pretrim_fastqc.zip") into pretrimFastQCResults

    """
    EXT="${fq}"
    EXT="\${EXT#*.}"
    ln -s "${fq}" "${fq.simpleName}_pretrim.\${EXT}"

    fastqc "${fq.simpleName}_pretrim.\${EXT}"
    """
}


pretrimFastQCResults.into {
    pretrimFastQC4SplitMultiQC;
    pretrimFastQC4JointMultiQC;
}


process cutadapt {
    label "python3"
    label "cutadapt"

    tag { bName }
    publishDir "${params.outdir}/trimmed"

    input:
    set val(bName), val(fwd_out), val(rev_out), file("fwd"), file("rev") from fastqPairs4Cutadapt
        .map { bName, fwd, rev -> [bName, fwd.name, rev.name, fwd, rev] }

    output:
    set val(bName), file(fwd_out), file(rev_out), file("*.txt") into trimmedPairs

    """
    # This bit is just so that we can keep the extensions for cutadapt
    # to autodetect the format.
    FWD="init_${fwd_out}"
    ln -sf fwd "\${FWD}"

    REV="init_${rev_out}"
    ln -sf rev "\${REV}"

    # Run cutadapt 3 times to trim that stuff yo.
    cutadapt \
        --quality-cutoff 25,25 \
        -a "${params.forward_adapter}" \
        -A "${params.reverse_adapter}" \
        --minimum-length 50 \
        -n 1 \
        --cores ${task.cpus} \
        -o "tmp1_${fwd_out}" \
        -p "tmp1_${rev_out}" \
        "\${FWD}" \
        "\${REV}" \
        > ${bName}_cutadapt_pass1.txt

    cutadapt \
        --quality-cutoff 25,25 \
        -a "${params.forward_adapter}" \
        -A "${params.reverse_adapter}" \
        --minimum-length 50 \
        -n 1 \
        --cores ${task.cpus} \
        -o "tmp2_${fwd_out}" \
        -p "tmp2_${rev_out}" \
        "tmp1_${fwd_out}" \
        "tmp1_${rev_out}" \
        > ${bName}_cutadapt_pass2.txt

    cutadapt \
        --quality-cutoff 25,25 \
        -a "${params.forward_adapter}" \
        -A "${params.reverse_adapter}" \
        --minimum-length 50 \
        -n 1 \
        --cores ${task.cpus} \
        -o "${fwd_out}" \
        -p "${rev_out}" \
        "tmp2_${fwd_out}" \
        "tmp2_${rev_out}" \
        > ${bName}_cutadapt_pass3.txt
    """
}


trimmedPairs.into {
    trimmedPairs4FastQC;
    trimmedPairs4MultiQC;
    trimmedPairs4UniqueKmers;
}


process posttrimFastQC {
    label "java"
    label "fastqc"
    tag { fq.simpleName }

    publishDir "${params.outdir}/read_stats"

    input:
    set val(bName), val(direction), file(fq) from trimmedPairs4FastQC
        .flatMap {bName, fwd, rev, txt -> [[bName, "forward", fwd], [bName, "reverse", rev]]}

    output:
    set val(bName), val(direction), file("${fq.simpleName}_posttrim_fastqc.html"),
        file("${fq.simpleName}_posttrim_fastqc.zip") into posttrimFastQCResults

    """
    EXT="${fq}"
    EXT="\${EXT#*.}"
    ln -s "${fq}" "${fq.simpleName}_posttrim.\${EXT}"

    fastqc "${fq.simpleName}_posttrim.\${EXT}"
    """
}


posttrimFastQCResults.into {
    posttrimFastQC4SplitMultiQC;
    posttrimFastQC4JointMultiQC;
}


process splitMultiQC {
    label "python3"
    label "multiqc"

    publishDir "${params.outdir}/read_stats"

    input:
    set val(direction), file("*"), file("*") from pretrimFastQC4SplitMultiQC
        .concat(posttrimFastQC4SplitMultiQC)
        .map {b, direction, html, zip -> [direction, html, zip]}
        .groupTuple(by: 0)

    output:
    set file("multiqc_${direction}.html"), file("multiqc_${direction}_data") into splitMultiQCResults

    """
    multiqc . --filename "multiqc_${direction}"
    """
}


process jointMultiQC {
    label "python3"
    label "multiqc"

    publishDir "${params.outdir}/read_stats"

    input:
    file "*" from pretrimFastQC4JointMultiQC
        .concat(posttrimFastQC4JointMultiQC)
        .flatMap {b, direction, html, zip -> [html, zip]}
        .collect()
    file "*" from trimmedPairs4MultiQC
        .flatMap {bName, fwd, rev, txt -> txt}
        .collect()

    output:
    set file("multiqc.html"), file("multiqc_data") into jointMultiQCResults

    """
    multiqc . --filename "multiqc"
    """
}

process uniqueKmers {
    label "python3"
    label "khmer"
    publishDir "${params.outdir}/read_stats"

    input:
    file "*" from trimmedPairs4UniqueKmers.map { b, f, r, t -> [f, r] }.flatten().collect()

    output:
    set file("unique_kmers_report.txt"), file("unique_kmers.txt") into uniqueKmersResults

    """
    FILES="\$(ls *)"
    unique-kmers.py -R unique_kmers_report.txt.tmp -e 0.01 -k 31 * 2> unique_kmers.txt

    echo "\${FILES}" > unique_kmers_report.txt
    cat unique_kmers_report.txt.tmp >> unique_kmers_report.txt

    rm unique_kmers_report.txt.tmp
    """
}

if ( params.references ) {

    process bowtieIndex {
        label "bowtie2"

        tag { reference.baseName }

        input:
        file reference from references

        output:
        file "${reference.baseName}" into referenceIndexes

        """
        mkdir -p "${reference.baseName}/index"
        bowtie2-build --threads ${task.cpus} -f "${reference}" "${reference.baseName}/index"
        """
    }

    process bowtieAlign {
        label "bowtie2"
        publishDir "${params.outdir}/aligned/${reference.name}"

        tag { "${reference}-${bName}" }

        input:
        set val(bName), file(fwd_read), file(rev_read) from fastqPairs4Bowtie
        each file(reference) from referenceIndexes

        output:
        set val(reference.name), val(bName), file("${bName}.sam"), file("${bName}_bowtie_stats.txt") into alignedReads

        """
        bowtie2 \
            --threads ${task.cpus} \
            --fr \
            -q \
            --phred33 \
            --fast-local \
            --dovetail \
            --minins 150 \
            --maxins 800 \
            -x "${reference}/index" \
            -1 "${fwd_read}" \
            -2 "${rev_read}" \
            -S "${bName}.sam" \
            2> "${bName}_bowtie_stats.txt"
        """
    }

    alignedReads.into {
        alignedReads4Stats;
        alignedReads4FastQC;
        alignedReads4MultiQC;
    }

    process alignmentStats {
        label "samtools"
        publishDir "${params.outdir}/aligned_stats/${reference}"

        tag { "${reference}-${bName}" }

        input:
        set val(reference), val(bName), file(sam), file(stats) from alignedReads4Stats

        output:
        set val(reference), val(bName), file("${bName}.idxstats"), file("${bName}.flagstat"), file("${bName}.stats") into samtoolsStats

        """
        samtools sort -O bam -o "${sam.baseName}.bam" "${sam}"
        samtools idxstats "${sam.baseName}.bam" > "${bName}.idxstats"
        samtools flagstat "${sam.baseName}.bam" > "${bName}.flagstat"
        samtools stats "${sam.baseName}.bam" > "${bName}.stats"
        """
    }

    process alignmentFastQC {
        label "fastqc"
        publishDir "${params.outdir}/aligned_stats/${reference}"

        tag { "${reference}-${bName}" }

        input:
        set val(reference), val(bName), file(sam), file(stats) from alignedReads4FastQC

        output:
        set val(reference), val(bName), file("${sam.simpleName}_fastqc.html"),
            file("${sam.simpleName}_fastqc.zip") into alignmentFastQCResults

        """
        fastqc ${sam}
        """
    }

    files4AlignmentMultiQC = alignmentFastQCResults
        .map {ref, bname, html, zip -> [ref, html, zip]}
        .groupTuple( by: 0 )
        .join(
            alignedReads4MultiQC
                .map {ref, bName, sam, stats -> [ref, stats]}
                .groupTuple( by: 0 ),
            remainder: true
        )
        .join(
            samtoolsStats
                .map {ref, bName, idxstats, flagstat, stats ->
                      [ref, idxstats, flagstat, stats]}
                .groupTuple( by: 0 ),
            remainder: true
        )

    process alignmentMultiQC {
        label "python3"
        label "multiqc"

        publishDir "${params.outdir}/aligned_stats/${ref}"

        input:
        // ref mqc_html   mqc_zip    bt2_stats  idxstats   flagstat   stats
        set val(ref), file("*"), file("*"), file("*"),
            file("*"), file("*"), file("*") from files4AlignmentMultiQC

        output:
        set file("multiqc.html"), file("multiqc_data") into alignmentMultiQCResults

        """
        multiqc . --filename "multiqc"
        """
    }
}

process searchKraken {
    label "kraken"

    publishDir "contaminants"

    tag { bName }

    input:
    file "krakendb" from krakenDB
    set val(bName), file(fwd_read), file(rev_read) from fastqPairs4Kraken

    output:
    set file("${bName}.tsv"), file("${bName}_report.txt") into krakenResults

    """
    kraken2 \
        --confidence 0.2 \
        --minimum-base-quality 25 \
        --paired \
        --output "${bName}.tsv" \
        --report "${bName}_report.txt" \
        --db krakendb \
        "${fwd_read}" \
        "${rev_read}"
    """
}
