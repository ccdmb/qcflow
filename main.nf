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

// Default parameter values

params.fastq = false
params.references = false
params.adapter1 = "data/adapters_truseq_fwd.fasta"
params.adapter2 = "data/adapters_truseq_rev.fasta"
params.synthetic_contaminants = "data/synthetic_contaminants.fasta"
params.quality_filter_phred = 25
params.quality_trim_phred = 10
params.minimum_read_length = 50
params.contaminants = false
params.nomap = false
params.nomerge = false
params.krakendb = false


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


if ( params.adapter1 && params.adapter2 ) {
    adapter1 = Channel.fromPath(params.adapter1, checkIfExists: true, type: "file")
        .collectFile(name: 'adapter1.fasta', newLine: true, sort: "deep")
        .first()

    adapter2 = Channel.fromPath(params.adapter2, checkIfExists: true, type: "file")
        .collectFile(name: 'adapter2.fasta', newLine: true, sort: "deep")
        .first()
} else {
    log.info "Hey I need some adapter sequences to trim please."
    log.info "Some TruSeq adapters are in the 'data' folder."
    exit 1
}


if ( params.synthetic_contaminants ) {
    syntheticContaminantSeqs = Channel.fromPath(params.synthetic_contaminants, checkIfExists: true, type: "file")
        .collectFile(name: "synthetic_contaminants.fasta", newLine: true, sort: "deep")
        .first()
} else {
    log.info "Hey I need some synthetic contaminant sequences to filter out please."
    log.info "I suggest using at least the 'data/synthetic_contaminants.fasta' file."
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


if ( params.krakendb ) {
    krakenDB = Channel.fromPath( params.krakendb, checkIfExists: true, type: "dir" ).first()
}


if ( params.contaminants ) {
    contaminants = Channel.fromPath(params.contaminants, checkIfExists: true, type: "file")
        .collectFile(name: "contaminants.fasta", newLine: true, sort: "deep")
        .first()
}


// END OF INPUT VALIDATION


fastqPairs.into {
    fastqPairs4RawStats;
    fastqPairs4QC;
    fastqPairs4AdapterTrimming;
    fastqPairs4Cutadapt;
    fastqPairs4Alignment;
    fastqPairs4Kraken;
}


/*
 * Initial QC stats to check if we improve things.
 */
process getRawStats {
    label "java"
    label "bbmap"
    label "medium_task"

    tag { base_name }

    publishDir "${params.outdir}/read_processing/raw"

    input:
    set val(base_name), file(fwd_read), file(rev_read) from fastqPairs4RawStats

    output:
    set val(base_name), file("*.txt") into rawStats

    """
    bbcountunique.sh \
      -Xmx${task.memory.toGiga()}g \
      in1="${fwd_read}" \
      in2="${rev_read}" \
      out="${base_name}_raw_count_unique.txt"

    kmercountmulti.sh \
      in1="${fwd_read}" \
      in2="${rev_read}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_raw_kmercountmulti.txt"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="${fwd_read}" \
      in2="${rev_read}" \
      stats="${base_name}_raw_stats.txt" \
      bhist="${base_name}_raw_bhist.txt" \
      qhist="${base_name}_raw_qhist.txt" \
      qchist="${base_name}_raw_qchist.txt" \
      aqhist="${base_name}_raw_aqhist.txt" \
      bqhist="${base_name}_raw_bqhist.txt" \
      lhist="${base_name}_raw_lhist.txt" \
      gchist="${base_name}_raw_gchist.txt" \
      gcbins="auto"
    """
}

rawStats
  .map { b, f -> [b, "raw", f] }
  .tap { rawStats4PerStepQC }
  .flatMap { b, s, f -> f }
  .set { rawStats4JointQC }


/*
 * Trim adapters using BBDuk.
 */
process adapterTrimming {
    label "java"
    label "bbmap"
    label "medium_task"

    tag { base_name }

    publishDir "${params.outdir}/read_processing/adapter_trimmed"

    input:
    file "fwd_adapters.fasta" from adapter1
    file "rev_adapters.fasta" from adapter2
    set val(base_name), val(fwd_name), val(rev_name), file("fwd"), file("rev") from fastqPairs4AdapterTrimming.map { b, f, r -> [b, f.name, r.name, f, r] }

    output:
    set val(base_name), val(fwd_name), val(rev_name), file(fwd_name), file(rev_name) into adapterTrimmed
    set val(base_name), file("*.txt") into adapterTrimmedStats

    """
    # This keeps the extensions around so we can autodetect the format.
    FWD="in_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="in_${rev_name}"
    ln -sf rev "\${REV}"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="\${FWD}" \
      in2="\${REV}" \
      out1="${fwd_name}" \
      out2="${rev_name}" \
      ref="fwd_adapters.fasta,rev_adapters.fasta" \
      stats="${base_name}_adapter_trimmed_stats.txt" \
      bhist="${base_name}_adapter_trimmed_bhist.txt" \
      qhist="${base_name}_adapter_trimmed_qhist.txt" \
      qchist="${base_name}_adapter_trimmed_qchist.txt" \
      aqhist="${base_name}_adapter_trimmed_aqhist.txt" \
      bqhist="${base_name}_adapter_trimmed_bqhist.txt" \
      lhist="${base_name}_adapter_trimmed_lhist.txt" \
      gchist="${base_name}_adapter_trimmed_gchist.txt" \
      gcbins="auto" \
      ktrim=r \
      k=23 \
      mink=11 \
      hdist=1 \
      minlength="${params.minimum_read_length}"

    kmercountmulti.sh \
      in1="${fwd_name}" \
      in2="${rev_name}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_adapter_trimmed_kmercountmulti.txt"
    """
}

adapterTrimmed.into {
    adapterTrimmed4QC;
    adapterTrimmed4QualityTrimming;
}

adapterTrimmedStats
  .map { b, f -> [b, "adapter_trimmed", f] }
  .tap { adapterTrimmedStats4PerStepQC }
  .flatMap { b, s, f -> f }
  .set { adapterTrimmedStats4JointQC }


/*
 * Trim adapter trimmed reads from BBDuk.
 */
process qualityTrimming {
    label "java"
    label "bbmap"
    label "medium_task"

    tag { base_name }

    publishDir "${params.outdir}/read_processing/quality_trimmed"

    input:
    set val(base_name),
        val(fwd_name), val(rev_name),
        file("fwd"), file("rev") from adapterTrimmed4QualityTrimming

    output:
    set val(base_name), val(fwd_name), val(rev_name),
        file(fwd_name), file(rev_name) into qualityTrimmed
    set val(base_name), file("*.txt") into qualityTrimmedStats

    """
    # This keeps the extensions around so we can autodetect the format.
    FWD="in_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="in_${rev_name}"
    ln -sf rev "\${REV}"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="\${FWD}" \
      in2="\${REV}" \
      out1="${fwd_name}" \
      out2="${rev_name}" \
      stats="${base_name}_quality_trimmed_stats.txt" \
      bhist="${base_name}_quality_trimmed_bhist.txt" \
      qhist="${base_name}_quality_trimmed_qhist.txt" \
      qchist="${base_name}_quality_trimmed_qchist.txt" \
      aqhist="${base_name}_quality_trimmed_aqhist.txt" \
      bqhist="${base_name}_quality_trimmed_bqhist.txt" \
      lhist="${base_name}_quality_trimmed_lhist.txt" \
      gchist="${base_name}_quality_trimmed_gchist.txt" \
      gcbins="auto" \
      qtrim=r \
      trimq="${params.quality_trim_phred}" \
      minlength="${params.minimum_read_length}"

    kmercountmulti.sh \
      in1="${fwd_name}" \
      in2="${rev_name}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_quality_trimmed_kmercountmulti.txt"
    """
}

qualityTrimmed.into {
    qualityTrimmed4SyntheticContaminantFilter;
    qualityTrimmed4QC;
}

qualityTrimmedStats
  .map { b, f -> [b, "quality_trimmed", f] }
  .tap { qualityTrimmedStats4PerStepQC }
  .flatMap { b, s, f -> f }
  .set { qualityTrimmedStats4JointQC }


/*
 * Combo adapter and quality trimming with cutadapt.
 */
process cutadapt {
    label "python3"
    label "cutadapt"
    label "medium_task"

    tag { base_name }

    publishDir "${params.outdir}/read_processing/cutadapt_trimmed"

    input:
    file "fwd_adapters.fasta" from adapter1
    file "rev_adapters.fasta" from adapter2
    set val(base_name), val(fwd_name), val(rev_name), file("fwd"), file("rev") from fastqPairs4Cutadapt.map { b, f, r -> [b, f.name, r.name, f, r] }

    output:
    set val(base_name), val(fwd_name), val(rev_name), file(fwd_name), file(rev_name) into cutadaptTrimmed
    set val(base_name), file("*.txt") into cutadaptStats

    """
    # This bit is just so that we can keep the extensions for cutadapt
    # to autodetect the format.
    FWD="init_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="init_${rev_name}"
    ln -sf rev "\${REV}"

    # Run cutadapt 3 times to trim that stuff yo.
    cutadapt \
        --quality-cutoff "${params.quality_trim_phred},${params.quality_trim_phred}" \
        -a "file:fwd_adapters.fasta" \
        -A "file:rev_adapters.fasta" \
        --minimum-length "${params.minimum_read_length}" \
        -n 3 \
        --cores ${task.cpus} \
        -o "tmp1_${fwd_name}" \
        -p "tmp1_${rev_name}" \
        "\${FWD}" \
        "\${REV}" \
        > ${base_name}_cutadapt_pass1.txt

    cutadapt \
        --quality-cutoff "${params.quality_trim_phred},${params.quality_trim_phred}" \
        -a "file:fwd_adapters.fasta" \
        -A "file:rev_adapters.fasta" \
        --minimum-length "${params.minimum_read_length}" \
        -n 3 \
        --cores ${task.cpus} \
        -o "${fwd_name}" \
        -p "${rev_name}" \
        "tmp1_${fwd_name}" \
        "tmp1_${rev_name}" \
        > ${base_name}_cutadapt_pass2.txt
    """
}

cutadaptTrimmed.into {
    cutadaptTrimmed4BBQC;
    cutadaptTrimmed4QC;
    cutadaptTrimmed4SyntheticContaminantFilter;
}

cutadaptStats
  .map { b, f -> [b, "cutadapt_trimmed", f] }
  .tap { cutadaptStats4PerStepQC }
  .flatMap { b, s, f -> f }
  .set { cutadaptStats4JointQC }


/*
 * Get BBDuk stats for Cutadapt trimmed reads.
 */
process getCutadaptBBQCStats {
    label "java"
    label "bbmap"
    label "medium_task"

    tag { base_name }

    publishDir "${params.outdir}/read_processing/cutadapt_trimmed"

    input:
    set val(base_name), val(fwd_name), val(rev_name), file(fwd_read), file(rev_read) from cutadaptTrimmed4BBQC

    output:
    set val(base_name), file("*.txt") into cutadaptBBQCStats

    """
    bbcountunique.sh \
      -Xmx${task.memory.toGiga()}g \
      in1="${fwd_read}" \
      in2="${rev_read}" \
      out="${base_name}_cutadapt_trimmed_count_unique.txt"

    kmercountmulti.sh \
      in1="${fwd_read}" \
      in2="${rev_read}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_cutadapt_trimmed_kmercountmulti.txt"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="${fwd_read}" \
      in2="${rev_read}" \
      stats="${base_name}_cutadapt_trimmed_stats.txt" \
      bhist="${base_name}_cutadapt_trimmed_bhist.txt" \
      qhist="${base_name}_cutadapt_trimmed_qhist.txt" \
      qchist="${base_name}_cutadapt_trimmed_qchist.txt" \
      aqhist="${base_name}_cutadapt_trimmed_aqhist.txt" \
      bqhist="${base_name}_cutadapt_trimmed_bqhist.txt" \
      lhist="${base_name}_cutadapt_trimmed_lhist.txt" \
      gchist="${base_name}_cutadapt_trimmed_gchist.txt" \
      gcbins="auto"
    """
}

cutadaptBBQCStats
  .map { b, f -> [b, "raw", f] }
  .tap { cutadaptBBQCStats4PerStepQC }
  .flatMap { b, s, f -> f }
  .set { cutadaptBBQCStats4JointQC }


/*
 * Filter PHiX and other common contaminants from reads.
 */
process syntheticContaminantFilter {
    label "java"
    label "bbmap"
    label "medium_task"

    tag { base_name }

    publishDir "${params.outdir}/read_processing/synthetic_contaminant_filtered"

    input:
    file "synthetic_contaminants.fasta" from syntheticContaminantSeqs
    set val(base_name), val(fwd_name), val(rev_name),
        file("fwd"), file("rev") from cutadaptTrimmed4SyntheticContaminantFilter

    output:
    set val(base_name), val(fwd_name), val(rev_name),
        file(fwd_name), file(rev_name) into syntheticContaminantFiltered
    set val(base_name), file("*.txt") into syntheticContaminantStats

    """
    # This keeps the extensions around so we can autodetect the format.
    FWD="in_${fwd_name}"
    ln -sf fwd "\${FWD}"

    REV="in_${rev_name}"
    ln -sf rev "\${REV}"

    bbduk.sh \
      -Xmx${task.memory.toGiga()}g \
      t=${task.cpus} \
      in1="\${FWD}" \
      in2="\${REV}" \
      out1="${fwd_name}" \
      out2="${rev_name}" \
      ref="synthetic_contaminants.fasta" \
      stats="${base_name}_synthetic_contaminant_filtered_stats.txt" \
      bhist="${base_name}_synthetic_contaminant_filtered_bhist.txt" \
      qhist="${base_name}_synthetic_contaminant_filtered_qhist.txt" \
      qchist="${base_name}_synthetic_contaminant_filtered_qchist.txt" \
      aqhist="${base_name}_synthetic_contaminant_filtered_aqhist.txt" \
      bqhist="${base_name}_synthetic_contaminant_filtered_bqhist.txt" \
      lhist="${base_name}_synthetic_contaminant_filtered_lhist.txt" \
      gchist="${base_name}_synthetic_contaminant_filtered_gchist.txt" \
      gcbins="auto" \
      k=31 \
      hdist=1 \
      mcf=0.5 \
      minavgquality="${params.quality_filter_phred}" \
      minlength="${params.minimum_read_length}"

    kmercountmulti.sh \
      in1="${fwd_name}" \
      in2="${rev_name}" \
      sweep=25,31,37,45,55,67,81,91 \
      stdev \
      out="${base_name}_synthetic_contaminant_filtered_kmercountmulti.txt"
    """
}

syntheticContaminantFiltered.into {
    syntheticContaminantFiltered4ContaminantFilter;
    syntheticContaminantFiltered4QC;
    syntheticContaminantFiltered4MergePairs;
}

syntheticContaminantStats
  .map { b, f -> [b, "synthetic_contaminant_filtered", f] }
  .tap { syntheticContaminantStats4PerStepQC }
  .flatMap { b, s, f -> f }
  .set { syntheticContaminantStats4JointQC }


if ( params.references && params.contaminants ) {

    /*
     * Mask the contaminants using the references.
     * This just converts all shared kmers to N's so that we avoid filtering out
     * genuine genetic content.
     */
    process maskContaminants {
        label "java"
        label "bbmap"
        label "medium_task"

        input:
        file "contaminants.fasta" from contaminants
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
}


if ( params.contaminants ) {

    /*
     * Filter contaminants from the reads.
     * Use pretty strict settings to avoid false positives.
     */
    process contaminantFilter {
        label "java"
        label "bbmap"
        label "bigmem_task" 

        tag { base_name }

        publishDir "${params.outdir}/read_processing/contaminant_filtered"

        input:
        file "contaminants.fasta" from maskedContaminants
        set val(base_name), val(fwd_name), val(rev_name),
            file("fwd"), file("rev") from syntheticContaminantFiltered4ContaminantFilter

        output:
        set val(base_name), val(fwd_name), val(rev_name),
            file(fwd_name), file(rev_name) into contaminantFiltered
        set val(base_name), file("*.txt") into contaminantStats

        """
        # This keeps the extensions around so we can autodetect the format.
        FWD="in_${fwd_name}"
        ln -sf fwd "\${FWD}"

        REV="in_${rev_name}"
        ln -sf rev "\${REV}"

        bbduk.sh \
          -Xmx${task.memory.toGiga()}g \
          t=${task.cpus} \
          in1="\${FWD}" \
          in2="\${REV}" \
          out1="${fwd_name}" \
          out2="${rev_name}" \
          ref="contaminants.fasta" \
          stats="${base_name}_contaminant_filtered_stats.txt" \
          bhist="${base_name}_contaminant_filtered_bhist.txt" \
          qhist="${base_name}_contaminant_filtered_qhist.txt" \
          qchist="${base_name}_contaminant_filtered_qchist.txt" \
          aqhist="${base_name}_contaminant_filtered_aqhist.txt" \
          bqhist="${base_name}_contaminant_filtered_bqhist.txt" \
          lhist="${base_name}_contaminant_filtered_lhist.txt" \
          gchist="${base_name}_contaminant_filtered_gchist.txt" \
          gcbins="auto" \
          k=25 \
          hdist=0 \
          mcf=0.6

        kmercountmulti.sh \
          in1="${fwd_name}" \
          in2="${rev_name}" \
          sweep=25,31,37,45,55,67,81,91 \
          stdev \
          out="${base_name}_contaminant_filtered_kmercountmulti.txt"
        """
    }

    contaminantFiltered.into {
        contaminantFiltered4QC;
        contaminantFiltered4Align;
        contaminantFiltered4MergePairs;
    }

    contaminantStats
      .map { b, f -> [b, "contaminant_filtered", f] }
      .tap { contaminantStats4PerStepQC }
      .flatMap { b, s, f -> f }
      .set { contaminantStats4JointQC }
}


/*
 * Here we set up some channels to run qc with.
 * It's a bit cumbersome because of all of the mapping, and because
 * some steps are optional.
 */

tmpFiles4PerStepQC = rawStats4PerStepQC.concat(
    adapterTrimmedStats4PerStepQC,
    syntheticContaminantStats4PerStepQC,
    qualityTrimmedStats4PerStepQC,
    cutadaptStats4PerStepQC,
    cutadaptBBQCStats4PerStepQC
)

tmpFiles4JointQC = rawStats4JointQC.concat(
    adapterTrimmedStats4JointQC,
    syntheticContaminantStats4JointQC,
    qualityTrimmedStats4JointQC,
    cutadaptStats4JointQC,
    cutadaptBBQCStats4JointQC
)

tmpPairs4ReadQC = fastqPairs4QC.flatMap {
        b, f, r -> [[b, "raw", "forward", f],
                    [b, "raw", "reverse", r]]
    }.concat(
    adapterTrimmed4QC.flatMap {
        b, fn, rn, ff, rf -> [[b, "adapter_trimmed", "forward", ff],
                              [b, "adapter_trimmed", "reverse", rf]]
    },
    syntheticContaminantFiltered4QC.flatMap {
        b, fn, rn, ff, rf -> [[b, "synthetic_contaminant_filtered", "forward", ff],
                              [b, "synthetic_contaminant_filtered", "reverse", rf]]
    },
    cutadaptTrimmed4QC.flatMap {
        b, fn, rn, ff, rf -> [[b, "cutadapt_trimmed", "forward", ff],
                              [b, "cutadapt_trimmed", "reverse", rf]]
    },
    qualityTrimmed4QC.flatMap {
        b, fn, rn, ff, rf -> [[b, "quality_trimmed", "forward", ff],
                              [b, "quality_trimmed", "reverse", rf]]
    }
)

/*
 * Because contaminant filtering is optional, we have to handle both cases.
 */
if ( params.contaminants ) {
    files4PerStepQC = tmpFiles4PerStepQC.concat(contaminantStats4PerStepQC)
    files4JointQC = tmpFiles4JointQC.concat(contaminantStats4JointQC)

    pairs4ReadQC = tmpPairs4ReadQC.concat(
        contaminantFiltered4QC.flatMap {
            b, fn, rn, ff, rf -> [
                [b, "contaminant_filtered", "forward", ff],
                [b, "contaminant_filtered", "reverse", rf]
            ]
        }
    )
} else {
    files4PerStepQC = tmpFiles4PerStepQC
    files4JointQC = tmpFiles4JointQC
    pairs4ReadQC = tmpPairs4ReadQC
}


/*
 * This just runs fastqc on each fastq file from each step individually.
 */
process readQC {
    label "java"
    label "fastqc"
    label "small_task"

    tag { "${base_name} - ${step} - ${direction}" }

    publishDir "${params.outdir}/read_processing/${step}"

    input:
    set val(base_name), val(step), val(direction), file(fq) from pairs4ReadQC

    output:
    set val(base_name), val(step), val(direction), file("${fq.simpleName}_${step}_fastqc.html"),
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


/*
 * Get multiqc files for each QC processing step individually.
 * Don't bother with fastqc because bbduks stats cover it.
 */
process perStepQC {
    label "python3"
    label "multiqc"
    label "small_task"

    tag { step }

    publishDir "${params.outdir}/read_processing/${step}"

    input:
    set val(step), file("*") from files4PerStepQC
        .flatMap { b, s, f -> f.collect {[s, it]} }
        .groupTuple(by: 0)

    output:
    set val(step), file("multiqc.html"), file("multiqc_data") into perStepQCResults

    """
    multiqc . --filename "multiqc"
    """
}


/*
 * Because fastqc processes each read separately and there's no way to
 * tell multiqc to plot R1 and R2 separately, we do two multiqc reports.
 */
process splitQC {
    label "python3"
    label "multiqc"
    label "small_task"

    tag { direction }

    publishDir "${params.outdir}/read_processing"

    input:
    set val(direction), file("*"), file("*") from readQC4SplitQC
        .map { b, step, direction, html, zip -> [direction, html, zip] }
        .groupTuple(by: 0)

    output:
    set file("multiqc_${direction}.html"), file("multiqc_${direction}_data") into splitMultiQCResults

    """
    multiqc . --filename "multiqc_${direction}"
    """
}


/*
 * Throw every QC result from every step and read orientation into one
 * MultiQC report.
 */
process jointQC {
    label "python3"
    label "multiqc"
    label "small_task"

    publishDir "${params.outdir}/read_processing"

    input:
    file "*" from readQC4JointQC
        .flatMap { b, step, direction, html, zip -> [html, zip] }
        .collect()
    file "*" from files4JointQC.collect()

    output:
    set file("multiqc.html"), file("multiqc_data") into jointMultiQCResults

    """
    multiqc . --filename "multiqc"
    """
}


if ( params.references && !params.nomap ) {

    /*
     * Construct genome indices for each reference to align to.
     */
    process genomeIndex {
        label "java"
        label "bbmap"
        label "biggish_task"

        tag { reference.baseName }

        input:
        file reference from references4Alignment

        output:
        set val(reference.baseName), file("ref") into referenceIndexes

        """
        bbmap.sh ref="${reference}"
        """
    }


    /*
     * Align reads to each reference and collect stats.
     */
    process genomeAlign {
        label "java"
        label "bbmap"
        label "biggish_task"

        publishDir "${params.outdir}/aligned/${ref_name}"

        tag { "${ref_name}-${base_name}" }

        input:
        set val(ref_name), file("ref"),
            val(base_name), file(fwd_read), file(rev_read) from referenceIndexes
            .combine(fastqPairs4Alignment)

        output:
        set val(ref_name), val(base_name), file("${base_name}.sam") into alignedReads
        set val(ref_name), file("*.txt") into alignedStats

        """
        bbmap.sh \
          -Xmx${task.memory.toGiga()}g \
          threads=${task.cpus} \
          in1="${fwd_read}" \
          in2="${rev_read}" \
          out="${base_name}.sam" \
          fast \
          local \
          covstats="${base_name}_constats.txt" \
          covhist="${base_name}_covhist.txt" \
          basecov="${base_name}_basecov.txt" \
          bincov="${base_name}_bincov.txt" \
          bhist="${base_name}_bhist.txt" \
          qhist="${base_name}_qhist.txt" \
          aqhist="${base_name}_aqhist.txt" \
          lhist="${base_name}_lhist.txt" \
          ihist="${base_name}_ihist.txt" \
          ehist="${base_name}_ehist.txt" \
          qahist="${base_name}_qahist.txt" \
          indelhist="${base_name}_indelhist.txt" \
          mhist="${base_name}_mhist.txt" \
          gchist="${base_name}_gchist.txt" \
          idhist="${base_name}_idhist.txt" \
          scafstats="${base_name}_scafstats.txt" \
          gcbins=auto \
          idbins=auto
        """
    }

    alignedReads.into {
        alignedReads4Stats;
        alignedReads4MultiQC;
    }


    /*
     * Get additional alignment statistics with samtools.
     */
    process alignmentStats {
        label "samtools"
        label "small_task"

        publishDir "${params.outdir}/aligned/${reference}"

        tag { "${reference}-${base_name}" }

        input:
        set val(reference), val(base_name), file(sam) from alignedReads4Stats

        output:
        set val(reference), file("${base_name}.idxstats"), file("${base_name}.flagstat"),
            file("${base_name}.stats") into samtoolsStats

        """
        samtools sort -O bam -o "${sam.baseName}.bam" "${sam}"
        samtools index "${sam.baseName}.bam"
        samtools idxstats "${sam.baseName}.bam" > "${base_name}.idxstats"
        samtools flagstat "${sam.baseName}.bam" > "${base_name}.flagstat"
        samtools stats "${sam.baseName}.bam" > "${base_name}.stats"

        rm -f ${sam.baseName}.{bam,bai}
        """
    }


    joined4AlignmentMultiQC = samtoolsStats
	.flatMap { r, i, f, s -> [[r, i], [r, f], [r, s]] }
	.concat(alignedStats.flatMap { r, s -> s.collect { [r, it] } })
	.filter { r, f -> !f.name.endsWith("qhist.txt") }
	.filter { r, f -> !f.name.endsWith("qahist.txt") }

    /*
     * Produce a multiqc report per reference for the isolates.
     */
    process alignmentMultiQC {
        label "python3"
        label "multiqc"
        label "small_task"

        publishDir "${params.outdir}/aligned/${ref}"

	tag { ref }

        input:
        set val(ref), file("*") from joined4AlignmentMultiQC
            .groupTuple(by: 0)

        output:
        set file("multiqc.html"), file("multiqc_data") into alignmentMultiQCResults

        """
        multiqc . --filename "multiqc"
        """
    }
}


if ( params.krakendb ) {

    /*
     * Classify the reads using Kraken to detect uncommon contamination.
     */
    process searchKraken {
        label "kraken"
        label "biggish_task"

        publishDir "contaminants"

        tag { base_name }

        input:
        file "krakendb" from krakenDB
        set val(base_name), file(fwd_read), file(rev_read) from fastqPairs4Kraken

        output:
        set file("${base_name}.tsv"), file("${base_name}_report.txt") into krakenResults

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


/*
 * Merge the fastq read channels and the empirical adapters.
 */

if ( params.contaminants ) {
    joined4MergePairs = contaminantFiltered4MergePairs
        .filter { b, t, fn, rn, ff, rf -> t == "untrimmed" }
        .map { b, t, fn, rn, ff, rf -> [b, fn, rn, ff, rf] }
} else {
    joined4MergePairs = syntheticContaminantFiltered4MergePairs
}


if ( !params.nomerge ) {
    /*
     * Brian Bushnell recommends loose for insert size estimation, and vstrict
     * for assembly.
     */
    stringencies = ['loose', 'vstrict']

    /*
     * Merge or stitch the read pairs.
     */
    process mergePairs {
        label "java"
        label "bbmap"
        label "big_task"

        publishDir "${params.outdir}/merged_reads/${stringency}"

        tag { "${base_name} - ${stringency}" }

        input:
        set val(base_name),
            val(fwd_name), val(rev_name),
            file("fwd"), file("rev") from joined4MergePairs
        each stringency from stringencies

        output:
        set val(base_name), file("${base_name}.fastq.gz"),
            file(fwd_name), file(rev_name) into mergedPairs
        set val(base_name), file("${base_name}_insert.txt") into mergedPairsInsert

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
          ${stringency} \
          rem \
          k=62 \
          extend2=50 \
          ecct
        """
    }
}
