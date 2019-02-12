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
      --adapters
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
params.adapters = false
params.synthetic_contaminants = false
params.no_quality_trim = false
params.quality_filter_phred = 25
params.quality_trim_phred = 20
params.contaminants = false
params.no_map = false
params.no_merge = false
params.nokraken = false
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


if ( params.adapters ) {
    adapterSeqs = Channel.fromPath(
        params.adapters,
        checkIfExists: true,
        type: "file"
    )
} else {
    log.info "Hey I need some adapter sequences to trim please."
    log.info "I suggest using 'adapters.fa' from bbmap, in the resources folder."
    exit 1
}

if ( params.synthetic_contaminants ) {
    syntheticContaminantSeqs = Channel.fromPath(
        params.synthetic_contaminants,
        checkIfExists: true,
        type: "file"
    )
} else {
    log.info "Hey I need some synthetic contaminant sequences to filter out please."
    log.info "I suggest using at least 'phix.fa' from bbmap, in the resources folder."
    exit 1
}

if ( params.contaminants ) {
    contaminants = Channel.fromPath(
        params.contaminants,
        checkIfExists: true,
        type: "file"
    )
} else {
    log.info "Hey I need some contaminant sequences to filter out please."
    log.info "I suggest using the refseq bacterial, viral, and a human genome."
    exit 1
}

if ( params.references ) {
    references = Channel.fromPath(
        params.references,
        checkIfExists: true,
        type: "file"
    )

    references.into {
        references4MaskContaminants;
        references4Alignment;
    }
}

if ( params.references && params.contaminants ) {
    process maskContaminants {
        label "java"
        label "bbmap"

        input:
        file "contaminants.fasta" from contaminants
            .collectFile(name: "contaminants.fasta", newLine: true, sort: "deep")
        file "references.fasta" from references4MaskContaminants
            .collectFile(name: "references.fasta", newLine: true, sort: "deep")

        output:
        file "masked_contaminants.fasta" into maskedContaminants

        """
        bbduk.sh \
          -Xmx${task.memory.toGiga()}g \
          t=${task.cpus} \
          in=contaminants.fasta \
          ref=references.fasta \
          out=masked_contaminants.fasta \
          k=25 \
          ktrim=N
        """
    }
} else if ( params.contaminants ) {
    maskedContaminants = contaminants
        .collectFile(name: "contaminants.fasta", newLine: true, sort: deep)
} else {
    log.error "Hey I reached a point in the code that I shouldn't be able to."
    log.error "I don't have contaminants or a reference."
    log.error "Please raise an issue."
    exit 1
}


if ( !params.krakendb && !params.nokraken) {

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

} else if ( !params.nokraken && file(params.krakendb).exists() ) {
    krakenDB = Channel.fromPath( params.krakendb )
} else if ( !params.nokraken ) {
    log.info "You specified a Kraken database, but it doesn't exist."
    log.info "Either provide the kraken database or disable kraken search with `--nokraken`."
    exit 1
}


fastqPairs.into {
    fastqPairs4QC;
    fastqPairs4SeparateNames;
    fastqPairs4FindAdapters;
    fastqPairs4Bowtie;
    fastqPairs4Kraken;
    fastqPairs4SplitBBDuk;
}

fastqPairs4SeparateNames
    .map { b, f, r -> [b, f.name, r.name, f, r]}
    .into {
        fastqPairs4AdapterTrimming;
        fastqPairs4Cutadapt;
    }

process findAdapters {
    label "java"
    label "bbmap"

    tag { base_name }

    input:
    set val(base_name), file(fwd_read), file(rev_read) from fastqPairs4FindAdapters

    output:
    set val(base_name), file("${fwd_read.baseName}_adapters.fasta"),
        file("${rev_read.baseName}_adapters.fasta") into foundAdapters

    """
    bbmerge.sh in="${fwd_read}" outa="${fwd_read.baseName}_adapters.fasta"
    bbmerge.sh in="${rev_read}" outa="${rev_read.baseName}_adapters.fasta"
    """
}

foundAdapters.into { foundAdapters4Trimming; foundAdapters4MergePairs }

joinedAdapters = fastqPairs4AdapterTrimming.join(foundAdapters4Trimming, by: 0)

process adapterTrimming {
    label "java"
    label "bbmap"

    tag { base_name }

    input:
    file "adapters.fasta" from adapterSeqs
        .collectFile(name: "adapters.fasta", newLine: true, sort: "deep")
    set val(base_name),
        val(fwd_name), val(rev_name), file("fwd"), file("rev"),
        file("fwd_adapters.fasta"), file("rev_adapters.fasta") from joinedAdapters

    output:
    set val(base_name),
        val(fwd_name), val(rev_name),
        file("${fwd_name}"), file("${rev_name}") into adapterTrimmed

    """
    # This keeps the extensions around so we can autodetect the format.
    FWD="in_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="in_${rev_name}"
    ln -sf rev "\${REV}"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="in_${fwd_name}" \
      in2="in_${rev_name}" \
      out1="${fwd_name}" \
      out2="${rev_name}" \
      ref="fwd_adapters.fasta,rev_adapters.fasta,adapters.fasta" \
      stats="${base_name}_stats.txt" \
      bhist="${base_name}_bhist.txt" \
      qhist="${base_name}_qhist.txt" \
      qchist="${base_name}_qchist.txt" \
      aqhist="${base_name}_aqhist.txt" \
      bqhist="${base_name}_bqhist.txt" \
      lhist="${base_name}_lhist.txt" \
      gchist="${base_name}_gchist.txt" \
      gcbins="auto" \
      ktrim=r \
      k=23 \
      mink=11 \
      hdist=1 \
      minlength=50 \
      minavgquality=25

    kmercountmulti.sh \
      in1="${fwd_name}" \
      in2="${rev_name}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_hist.txt"
    """
}

adapterTrimmed.into {
    adapterTrimmed4QC;
    adapterTrimmed4SyntheticContaminantFilter;
}

process syntheticContaminantFilter {
    label "java"
    label "bbmap"

    tag { base_name }

    input:
    file "synthetic_contaminants.fasta" from syntheticContaminantSeqs
        .collectFile(name: "synthetic_contaminants.fasta", newLine: true, sort: "deep")
    set val(base_name),
        val(fwd_name), val(rev_name),
        file("fwd"), file("rev") from adapterTrimmed4SyntheticContaminantFilter

    output:
    set val(base_name),
        val(fwd_name), val(rev_name),
        file("${fwd_name}"), file("${rev_name}") into syntheticContaminantFiltered
    set val(base_name), file("${base_name}_stats.txt") into syntheticContaminantStats

    """
    # This keeps the extensions around so we can autodetect the format.
    FWD="in_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="in_${rev_name}"
    ln -sf rev "\${REV}"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="in_${fwd_name}" \
      in2="in_${rev_name}" \
      out1="${fwd_name}" \
      out2="${rev_name}" \
      ref="synthetic_contaminants.fasta" \
      stats="${base_name}_stats.txt" \
      bhist="${base_name}_bhist.txt" \
      qhist="${base_name}_qhist.txt" \
      qchist="${base_name}_qchist.txt" \
      aqhist="${base_name}_aqhist.txt" \
      bqhist="${base_name}_bqhist.txt" \
      lhist="${base_name}_lhist.txt" \
      gchist="${base_name}_gchist.txt" \
      gcbins="auto" \
      k=31 \
      hdist=1 \
      mcf=0.4

    kmercountmulti.sh \
      in1="${fwd_name}" \
      in2="${rev_name}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_hist.txt"
    """
}

syntheticContaminantFiltered.into {
    syntheticContaminantFiltered4QualityTrimming;
    syntheticContaminantFiltered4ContaminantFilter;
    syntheticContaminantFiltered4QC;
}


if ( !params.no_quality_trim ) {
    process qualityTrimming {
        label "java"
        label "bbmap"

        tag { base_name }

        input:
        set val(base_name),
            val(fwd_name), val(rev_name),
            file("fwd"), file("rev") from syntheticContaminantFiltered4QualityTrimming

        output:
        set val(base_name),
            val(fwd_name), val(rev_name),
            file("${fwd_name}"), file("${rev_name}") into qualityTrimmed

        """
        # This keeps the extensions around so we can autodetect the format.
        FWD="in_${fwd_name}"
        ln -sf fwd "\${FWD}"

        REV="in_${rev_name}"
        ln -sf rev "\${REV}"

        bbduk.sh \
          -Xmx${task.memory.toGiga()}g \
          t=${task.cpus} \
          in1="in_${fwd_name}" \
          in2="in_${rev_name}" \
          out1="${fwd_name}" \
          out2="${rev_name}" \
          stats="${base_name}_stats.txt" \
          bhist="${base_name}_bhist.txt" \
          qhist="${base_name}_qhist.txt" \
          qchist="${base_name}_qchist.txt" \
          aqhist="${base_name}_aqhist.txt" \
          bqhist="${base_name}_bqhist.txt" \
          lhist="${base_name}_lhist.txt" \
          gchist="${base_name}_gchist.txt" \
          gcbins="auto" \
          maq="${params.quality_filter_phred}" \
          qtrim=r \
          trimq="${params.quality_trim_phred}" \
          minlength=50 \
          minavgquality=25

        kmercountmulti.sh \
          in1="${fwd_name}" \
          in2="${rev_name}" \
          sweep=25,31,37,45,55,67,81,91 \
          stdev \
          out="${base_name}_hist.txt"
        """
    }

    qualityTrimmed.into {
        qualityTrimmed4ContaminantFilter;
        qualityTrimmed4QC;
    }
}

// This just adds a new val field to differentiate trimmed from untrimmed output.
pairs4ContaminantFilter = syntheticContaminantFiltered4ContaminantFilter
    .map { b, fn, rn, ff, fr -> [b, "untrimmed", fn, rn, ff, fr] }
    .concat(
        qualityTrimmed4ContaminantFilter
            .map { b, fn, rn, ff, fr -> [b, "trimmed", fn, rn, ff, fr] }
    )

process contaminantFilter {
    label "java"
    label "bbmap"
    publishDir "${params.outdir}/bbduk_trimmed/${trimmed}"

    tag { "${base_name} - ${trimmed}" }

    input:
    file "contaminants.fasta" from maskedContaminants
    set val(base_name), val(trimmed),
        val(fwd_name), val(rev_name),
        file("fwd"), file("rev") from pairs4ContaminantFilter

    output:
    set val(base_name), val(trimmed),
        val(fwd_name), val(rev_name),
        file("${fwd_name}"), file("${rev_name}") into contaminantFiltered
    set val(base_name), file("${base_name}_stats.txt") into contaminantStats

    """
    # This keeps the extensions around so we can autodetect the format.
    FWD="in_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="in_${rev_name}"
    ln -sf rev "\${REV}"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="in_${fwd_name}" \
      in2="in_${rev_name}" \
      out1="${fwd_name}" \
      out2="${rev_name}" \
      ref="contaminants.fasta" \
      stats="${base_name}_stats.txt" \
      bhist="${base_name}_bhist.txt" \
      qhist="${base_name}_qhist.txt" \
      qchist="${base_name}_qchist.txt" \
      aqhist="${base_name}_aqhist.txt" \
      bqhist="${base_name}_bqhist.txt" \
      lhist="${base_name}_lhist.txt" \
      gchist="${base_name}_gchist.txt" \
      gcbins="auto" \
      k=27 \
      hdist=1 \
      mcf=0.6

    kmercountmulti.sh \
      in1="${fwd_name}" \
      in2="${rev_name}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_hist.txt"
    """
}

contaminantFiltered.into {
    contaminantFiltered4QC;
    contaminantFiltered4Align;
    contaminantFiltered4MergePairs;
}

process cutadapt {
    label "python3"
    label "cutadapt"

    tag { base_name }
    publishDir "${params.outdir}/cutadapt_trimmed"

    input:
    set val(base_name),
        val(fwd_out), val(rev_out),
        file("fwd"), file("rev") from fastqPairs4Cutadapt

    output:
    set val(base_name), file(fwd_out), file(rev_out), file("*.txt") into cutadaptTrimmed

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
        > ${base_name}_cutadapt_pass1.txt

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
        > ${base_name}_cutadapt_pass2.txt

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
        > ${base_name}_cutadapt_pass3.txt
    """
}

cutadaptTrimmed.into {
    cutadaptTrimmed4QC;
    cutadaptTrimmed4MultiQC;
}


pairs4ReadQC = Channel.create().concat(
    fastqPairs4QC.flatMap {
        b, f, r -> [[b, "raw", "forward", f],
                    [b, "raw", "reverse", r]]
    },
    adapterTrimmed4QC.flatMap {
        b, fn, rn, ff, rf -> [[b, "adapter_trimmed", "forward", ff],
                              [b, "adapter_trimmed", "reverse", rf]]
    },
    syntheticContaminantFiltered4QC.flatMap {
        b, fn, rn, ff, rf -> [[b, "synthetic_contaminant_filtered", "forward", ff],
                              [b, "synthetic_contaminant_filtered", "reverse", rf]]
    },
    cutadaptTrimmed4QC.flatMap {
        b, f, r, s -> [[b, "cutadapt_trimmed", "forward", f],
                       [b, "cutadapt_trimmed", "reverse", r]]
    },
    contaminantFiltered4QC.flatMap {
        b, t, fn, rn, ff, rf -> [[b, "contaminant_filtered_${t}", "forward", ff],
                                 [b, "contaminant_filtered_${t}", "forward", ff]]
    }
)

process readQC {
    label "java"
    label "fastqc"
    tag { "${base_name} - ${step} - ${direction}" }

    publishDir "${params.outdir}/read_stats"

    input:
    set val(base_name), val(step), val(direction), file(fq) from pairs4ReadQC

    output:
    set val(base_name), val(direction), file("${fq.simpleName}_${step}_fastqc.html"),
        file("${fq.simpleName}_${step}_fastqc.zip") into readQCResults

    """
    EXT="${fq}"
    EXT="\${EXT#*.}"
    ln -s "${fq}" "${fq.simpleName}_${step}.\${EXT}"

    fastqc "${fq.simpleName}_${step}.\${EXT}"
    """
}

readQCResults.into {
    readQC4SplitQC;
    readQC4JointQC;
}


process splitQC {
    label "python3"
    label "multiqc"

    publishDir "${params.outdir}/read_stats"

    input:
    set val(direction), file("*"), file("*") from readQC4SplitQC
        .map {b, direction, html, zip -> [direction, html, zip]}
        .groupTuple(by: 0)

    output:
    set file("multiqc_${direction}.html"), file("multiqc_${direction}_data") into splitMultiQCResults

    """
    multiqc . --filename "multiqc_${direction}"
    """
}


process jointQC {
    label "python3"
    label "multiqc"

    publishDir "${params.outdir}/read_stats"

    input:
    file "*" from readQC4JointQC
        .flatMap {b, direction, html, zip -> [html, zip]}
        .collect()
    file "*" from cutadaptTrimmed4MultiQC
        .flatMap { b, f, r, stats -> stats}
        .collect()

    output:
    set file("multiqc.html"), file("multiqc_data") into jointMultiQCResults

    """
    multiqc . --filename "multiqc"
    """
}


if ( params.references && !params.no_map ) {

    process bowtieIndex {
        label "bowtie2"

        tag { reference.baseName }

        input:
        file reference from references4Alignment

        output:
        file "${reference.baseName}" into referenceIndexes

        """
        mkdir -p "${reference.baseName}/index"
        bowtie2-build \
          --threads ${task.cpus} \
          -f "${reference}" \
          "${reference.baseName}/index"
        """
    }

    process bowtieAlign {
        label "bowtie2"
        publishDir "${params.outdir}/aligned/${reference.name}"

        tag { "${reference}-${base_name}" }

        input:
        set val(base_name), file(fwd_read), file(rev_read) from fastqPairs4Bowtie
        each file(reference) from referenceIndexes

        output:
        set val(reference.name), val(base_name), file("${base_name}.sam"),
            file("${base_name}_bowtie_stats.txt") into alignedReads

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
            -S "${base_name}.sam" \
            2> "${base_name}_bowtie_stats.txt"
        """
    }

    alignedReads.into {
        alignedReads4Stats;
        alignedReads4MultiQC;
    }

    process alignmentStats {
        label "samtools"
        publishDir "${params.outdir}/aligned_stats/${reference}"

        tag { "${reference}-${base_name}" }

        input:
        set val(reference), val(base_name), file(sam), file(stats) from alignedReads4Stats

        output:
        set val(reference), val(base_name), file("${base_name}.idxstats"),
            file("${base_name}.flagstat"), file("${base_name}.stats") into samtoolsStats

        """
        samtools sort -O bam -o "${sam.baseName}.bam" "${sam}"
        samtools index "${sam.baseName}.bam"
        samtools idxstats "${sam.baseName}.bam" > "${base_name}.idxstats"
        samtools flagstat "${sam.baseName}.bam" > "${base_name}.flagstat"
        samtools stats "${sam.baseName}.bam" > "${base_name}.stats"
        """
    }

    files4AlignmentMultiQC = alignedReads4MultiQC
            .map {ref, base_name, sam, stats -> [ref, stats]}
            .groupTuple( by: 0 )
        .join(
            samtoolsStats
                .map {ref, base_name, idxstats, flagstat, stats ->
                      [ref, idxstats, flagstat, stats]}
                .groupTuple( by: 0 ),
            remainder: true
        )

    process alignmentMultiQC {
        label "python3"
        label "multiqc"

        publishDir "${params.outdir}/aligned_stats/${ref}"

        input:
        // ref    bt2_stats  idxstats   flagstat   stats
        set val(ref), file("*"), file("*"), file("*"), file("*") from files4AlignmentMultiQC

        output:
        set file("multiqc.html"), file("multiqc_data") into alignmentMultiQCResults

        """
        multiqc . --filename "multiqc"
        """
    }
}


if ( !params.nokraken ) {
    process searchKraken {
        label "kraken"

        publishDir "contaminants"

        tag { base_name }

        input:
        file "krakendb" from krakenDB
        set val(base_name), file(fwd_read), file(rev_read) from fastqPairs4Kraken

        output:
        set file("${base_name}.tsv"), file("${base_name}_report.txt") into krakenResults

        when:
        !params.nokraken

        """
        kraken2 \
          --confidence 0.2 \
          --minimum-base-quality 25 \
          --paired \
          --output "${base_name}.tsv" \
          --report "${base_name}_report.txt" \
          --db krakendb \
          "${fwd_read}" \
          "${rev_read}"
        """
    }
}

joined4MergePairs = contaminantFiltered4MergePairs
    .filter { b, t, fn, rn, ff, rf -> t == "untrimmed" }
    .map { b, t, fn, rn, ff, rf -> [b, fn, rn, ff, rf] }
    .join(foundAdapters4MergePairs, by: 0)

stringencies = ['loose', 'vstrict']

process mergePairs {
    label "java"
    label "bbmap"

    publishDir "${params.outdir}/merged_reads/${stringency}"

    tag { "${base_name} - ${stringency}" }

    input:
    set val(base_name),
        val(fwd_name), val(rev_name),
        file("fwd"), file("rev"),
        file("fwd_adapters.fasta"), file("rev_adapters.fasta") from joined4MergePairs
    each stringency from stringencies

    output:
    set val(base_name), file("${base_name}.fastq.gz"),
        file(fwd_name), file(rev_name) into mergedPairs
    set val(base_name), file("${base_name}_insert.txt") into mergedPairsInsert

    when:
    !params.no_merge

    """
    # This keeps the extensions around so we can autodetect the format.
    FWD="in_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="in_${rev_name}"
    ln -sf rev "\${REV}"

    bbmerge-auto.sh \
      in1="\${FWD}" \
      in2="\${REV}" \
      out="${base_name}.fastq.gz" \
      outu1="${fwd_name}" \
      outu2="${rev_name}" \
      ihist=${base_name}_insert.txt \
      adapter1="fwd_adapters.fasta" \
      adapter2="rev_adapters.fasta" \
      ${stringency} \
      rem \
      k=62 \
      extend2=50 \
      ecct
    """
}
