#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # qcflow

    A [nextflow](www.nextflow.io) pipeline for running multiple illumina
    short-read filtering, trimming, and QC steps. QC is performed at each stage.

    ## Usage

    ```bash
    nextflow run -resume darcyabjones/qcflow \
      --fastq "fastq/*_R{1,2}.fastq.gz" \
      --adapter1 "data/adapters_truseq_fwd.fasta" \
      --adapter2 "data/adapters_truseq_rev.fasta" \
      --synthetic_contaminants "data/synthetic_contaminants.fasta" \
      --references "my_genome*.fasta" \
      --map
    ```

    Note that the quotes around globbing patterns are important for some
    shells that automatically expand globs (e.g. zsh).

    ## Mandatory Arguments

    ```
    param                   | description
    ---------------------------------------------------------------------------
    `--fastq <fastq pairs>` | The fastq pairs to process. This must be
                            | provided as a glob pattern to capture the pairs.
                            | E.G. `sample{1,2}.fastq` will capture
                            | `sample1.fastq sample2.fastq`.
    ```

    ## Options

    ```
    param                      | default          | description
    ---------------------------------------------------------------------------
    `--outdir <path>`          | ./               | The base directory to
                               |                  | publish the results in.

    `--references <*fasta>`    | none             | A glob pattern of reference
                               |                  | genomes to use for mapping
                               |                  | and for masking contaminant
                               |                  | filtering.

    `--adapter1 <*fasta>`      | data/            | Adapter sequences to trim
                               | truseq_fwd.fasta | from the 5' end.

    `--adapter2 <*fasta>`      | data/            | Adapter sequences to trim
                               | truseq_rev.fasta | from the 3' end.

    `--scontaminants <*fasta>` | data/            | Synthetic contaminants to
                               | synt_cont.fasta  | filter from the reads.
                               |                  | This is for things like
                               |                  | PHiX or primer dimer, or
                               |                  | common lab vectors that are
                               |                  | likely contaminants of
                               |                  | sequencing rather than of
                               |                  | the samples themselves.

    `--contaminants <*fasta>   | none             | Filter reads that match
                               |                  | these sequences out.
                               |                  | This is for filtering out
                               |                  | sample contaminants, e.g.
                               |                  | bacteria or endophytes.
                               |                  | If a reference genome is
                               |                  | provided, regions in this
                               |                  | database matching the
                               |                  | reference will be masked.
                               |                  | Generally I would run kraken
                               |                  | to check for contamination
                               |                  | before using this.

    `--filter_phred <int>`     | 5                | Filter out reads that have
                               |                  | lower average phred scores
                               |                  | than this. Keep this low.
                               |                  | bbduk seems to be a bit
                               |                  | "filter happy" for this.

    `--trim_phred <int>`       | 2                | Trim bases with phred
                               |                  | qualities lower than this
                               |                  | off the end of reads.
                               |                  | Generally quality trimming
                               |                  | is not useful unless you
                               |                  | have very poor data.

    `--use_bbduk_trim`         | false            | Use bbduk for trimming
                               |                  | instead of cutadapt.

    `--min_read_length <int>`  | 50               | Filter out reads with
                               |                  | lengths less than this
                               |                  | after trimming.

    `--map`                    | false            | Align the raw reads to the
                               |                  | reference genomes and get
                               |                  | qc stats. Useful for
                               |                  | estimating insert/fragment
                               |                  | size, or error rates.

    `--merge`                  | false            | Merge paired end reads using
                               |                  | bbmerge. This is intended
                               |                  | for fragment size estimation
                               |                  | when a reference is
                               |                  | unavailable. The parameters
                               |                  | are not appropriate for
                               |                  | merging prior to assembly.

    `--krakendb <dir>`         | none             | Search for matches to
                               |                  | potential contaminants in
                               |                  | this kraken database.
                               |                  | This should be seen as a
                               |                  | first pass, check. If you
                               |                  | see something more
                               |                  | substantial you might want
                               |                  | to do a more accurate
                               |                  | alignment.

    `--kraken_low_mem`         | false            | Prevents kraken from
                               |                  | memory mapping the file.
                               |                  | This is used in two contexts.
                               |                  | 1) preserve ram, run slow.
                               |                  | 2) See "running kraken".
    ```

    ## Running kraken

    Kraken is useful as a first pass to detect potential contaminants in your
    sequencing. Which you might then choose to filter out using something like
    bbduk or bbsplit. It should be noted that the taxonomic assignments given
    by kraken are highly dependent of the database you provide, so you should
    usually try aligning with a more sensitive aligner (e.g. bma-mem, bbmap),
    or database search tool (e.g. blast+) before pushing the contamination red
    button.

    There is a slurm script to download and prepare a kraken database in the
    `batch_scripts` directory in the repo.

    Kraken is by far the slowest step in the pipeline if it isn't set up well.
    The issue is that for every fastq pair, it tries to load the database into
    ram again, which is slow. The best thing to do is to copy the database to
    a memory mounted filesystem before running the pipelineand use the
    `--kraken_low_mem` option. This way the file stays in memory for the full
    time and kraken won't try to load it into memory again.
    Many super computers already have memory mounted drives on nodes.
    E.G. pawsey mounts `/tmp` in ram, so you can just copy it there as part
    of the job submission script.

    ## Outputs

    There's a lot. I'll describe it one day.

    ## Requirements

    * `BBMap` <https://sourceforge.net/projects/bbmap/>.
      Developed with v38.39.
    * `fastqc` <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>.
      Developed with v0.11.8.
    * `cutadapt` <https://cutadapt.readthedocs.io/en/stable/guide.html>.
      Developed with v1.18.
      Required by default, but can be made optional with `--use_bbduk_trim`.
    * `multiqc` <https://multiqc.info/>
      Developed with v1.7
    * `kraken2` <https://ccb.jhu.edu/software/kraken2/>.
      Developed with v2.0.7-beta.
      Optional, only required for Kraken steps.
    * `samtools` <http://www.htslib.org/doc/samtools.html>.
      Developed with v1.6.
      Optional, only needed if running the mapping steps.
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

// Default parameter values

params.fastq = false
params.references = false
params.adapter1 = "data/truseq_fwd.fasta"
params.adapter2 = "data/truseq_rev.fasta"
params.scontaminants = "data/synth_cont.fasta"
params.filter_phred = 5
params.trim_phred = 2
params.min_read_length = 50
params.contaminants = false
params.map = false
params.merge = false
params.krakendb = false
params.kraken_low_mem = false
params.use_bbduk_trim = false


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
    adapter1 = Channel
        .fromPath(params.adapter1, checkIfExists: true, type: "file")
        .collectFile(name: 'adapter1.fasta', newLine: true, sort: "deep")
        .first()

    adapter2 = Channel
        .fromPath(params.adapter2, checkIfExists: true, type: "file")
        .collectFile(name: 'adapter2.fasta', newLine: true, sort: "deep")
        .first()
} else {
    log.info "Hey I need some adapter sequences to trim please."
    log.info "Some TruSeq adapters are in the 'data' folder."
    exit 1
}


if ( params.scontaminants ) {
    syntheticContaminantSeqs = Channel
        .fromPath(params.scontaminants, checkIfExists: true, type: "file")
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
    krakenDB = Channel
        .fromPath( params.krakendb, checkIfExists: true, type: "dir" )
        .first()
}


if ( params.contaminants ) {
    contaminants = Channel
        .fromPath(params.contaminants, checkIfExists: true, type: "file")
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


if ( params.use_bbduk_trim ) {
    /*
     * Trim adapters using BBDuk.
     */
    process adapterTrimming {
        label "bbmap"
        label "medium_task"

        tag { base_name }

        publishDir "${params.outdir}/read_processing/adapter_trimmed"

        input:
        file "fwd_adapters.fasta" from adapter1
        file "rev_adapters.fasta" from adapter2
        set val(base_name), val(fwd_name), val(rev_name),
            file("fwd"), file("rev") from fastqPairs4AdapterTrimming
                .map { b, f, r -> [b, f.name, r.name, f, r] }

        output:
        set val(base_name), val(fwd_name), val(rev_name),
            file(fwd_name), file(rev_name) into adapterTrimmed
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
          minlength="${params.min_read_length}"

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
        label "bbmap"
        label "medium_task"

        tag { base_name }

        publishDir "${params.outdir}/read_processing/quality_trimmed"

        input:
        set val(base_name), val(fwd_name), val(rev_name),
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
          trimq="${params.trim_phred}" \
          minlength="${params.min_read_length}"

        kmercountmulti.sh \
          in1="${fwd_name}" \
          in2="${rev_name}" \
          sweep=25,31,37,45,55,67,81,91 \
          stdev \
          out="${base_name}_quality_trimmed_kmercountmulti.txt"
        """
    }

    qualityTrimmed.into {
        reads4SyntheticContaminantFilter;
        qualityTrimmed4QC;
    }

    qualityTrimmedStats
        .map { b, f -> [b, "quality_trimmed", f] }
        .tap { qualityTrimmedStats4PerStepQC }
        .flatMap { b, s, f -> f }
        .set { qualityTrimmedStats4JointQC }

} else {
    /*
     * Combo adapter and quality trimming with cutadapt.
     */
    process cutadapt {
        label "cutadapt"
        label "medium_task"

        tag { base_name }

        publishDir "${params.outdir}/read_processing/cutadapt_trimmed"

        input:
        file "fwd_adapters.fasta" from adapter1
        file "rev_adapters.fasta" from adapter2
        set val(base_name), val(fwd_name), val(rev_name),
            file("fwd"), file("rev") from fastqPairs4Cutadapt
                .map { b, f, r -> [b, f.name, r.name, f, r] }

        output:
        set val(base_name), val(fwd_name), val(rev_name),
            file(fwd_name), file(rev_name) into cutadaptTrimmed
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
            --quality-cutoff "${params.trim_phred},${params.trim_phred}" \
            -a "file:fwd_adapters.fasta" \
            -A "file:rev_adapters.fasta" \
            --minimum-length "${params.min_read_length}" \
            -n 3 \
            --cores ${task.cpus} \
            -o "tmp1_${fwd_name}" \
            -p "tmp1_${rev_name}" \
            "\${FWD}" \
            "\${REV}" \
            > ${base_name}_cutadapt_pass1.txt

        cutadapt \
            --quality-cutoff "${params.trim_phred},${params.trim_phred}" \
            -a "file:fwd_adapters.fasta" \
            -A "file:rev_adapters.fasta" \
            --minimum-length "${params.min_read_length}" \
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
        reads4SyntheticContaminantFilter;
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
        label "bbmap"
        label "medium_task"

        tag { base_name }

        publishDir "${params.outdir}/read_processing/cutadapt_trimmed"

        input:
        set val(base_name), val(fwd_name), val(rev_name),
            file(fwd_read), file(rev_read) from cutadaptTrimmed4BBQC

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

} // ENDIF params.use_bbduk_trim


/*
 * Filter PHiX and other common contaminants from reads.
 */
process syntheticContaminantFilter {
    label "bbmap"
    label "medium_task"

    tag { base_name }

    publishDir "${params.outdir}/read_processing/synthetic_contaminant_filtered"

    input:
    file "synthetic_contaminants.fasta" from syntheticContaminantSeqs
    set val(base_name), val(fwd_name), val(rev_name),
        file("fwd"), file("rev") from reads4SyntheticContaminantFilter

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
      hdist=0 \
      mcf=0.7 \
      minavgquality="${params.filter_phred}" \
      minlength="${params.min_read_length}"

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
          k=27 \
          hdist=0 \
          mcf=0.7

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
if ( params.use_bbduk_trim ) {
    tmpFiles4PerStepQC = rawStats4PerStepQC.concat(
        adapterTrimmedStats4PerStepQC,
        syntheticContaminantStats4PerStepQC,
        qualityTrimmedStats4PerStepQC,
    )

    tmpFiles4JointQC = rawStats4JointQC.concat(
        adapterTrimmedStats4JointQC,
        syntheticContaminantStats4JointQC,
        qualityTrimmedStats4JointQC,
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
        qualityTrimmed4QC.flatMap {
            b, fn, rn, ff, rf -> [[b, "quality_trimmed", "forward", ff],
                                  [b, "quality_trimmed", "reverse", rf]]
        }
    )

} else {
    tmpFiles4PerStepQC = rawStats4PerStepQC.concat(
        syntheticContaminantStats4PerStepQC,
        cutadaptStats4PerStepQC,
        cutadaptBBQCStats4PerStepQC
    )

    tmpFiles4JointQC = rawStats4JointQC.concat(
        syntheticContaminantStats4JointQC,
        cutadaptStats4JointQC,
        cutadaptBBQCStats4JointQC
    )

    tmpPairs4ReadQC = fastqPairs4QC.flatMap {
            b, f, r -> [[b, "raw", "forward", f],
                        [b, "raw", "reverse", r]]
        }.concat(
        syntheticContaminantFiltered4QC.flatMap {
            b, fn, rn, ff, rf -> [[b, "synthetic_contaminant_filtered", "forward", ff],
                                  [b, "synthetic_contaminant_filtered", "reverse", rf]]
        },
        cutadaptTrimmed4QC.flatMap {
            b, fn, rn, ff, rf -> [[b, "cutadapt_trimmed", "forward", ff],
                                  [b, "cutadapt_trimmed", "reverse", rf]]
        }
    )
}


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
    label "fastqc"
    label "small_task"

    tag { "${base_name} - ${step} - ${direction}" }

    publishDir "${params.outdir}/read_processing/${step}"

    input:
    set val(base_name), val(step), val(direction), file(fq) from pairs4ReadQC

    output:
    set val(base_name), val(step), val(direction),
        file("${fq.simpleName}_${step}_fastqc.html"),
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
    label "multiqc"
    label "small_task"

    tag { direction }

    publishDir "${params.outdir}/read_processing"

    input:
    set val(direction), file("*"), file("*") from readQC4SplitQC
        .map { b, step, direction, html, zip -> [direction, html, zip] }
        .groupTuple(by: 0)

    output:
    set file("multiqc_${direction}.html"),
        file("multiqc_${direction}_data") into splitMultiQCResults

    """
    multiqc . --filename "multiqc_${direction}"
    """
}


/*
 * Throw every QC result from every step and read orientation into one
 * MultiQC report.
 */
process jointQC {
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


if ( params.references && params.map ) {

    /*
     * Construct genome indices for each reference to align to.
     */
    process genomeIndex {
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
        label "big_task"

        publishDir "contaminants"

        tag { base_name }

        input:
        file "krakendb" from krakenDB
        set val(base_name), file(fwd_read), file(rev_read) from fastqPairs4Kraken

        output:
        set file("${base_name}.tsv"), file("${base_name}_report.txt") into krakenResults

    script:
    if ( params.kraken_low_mem ) {
        mmap = "--memory-mapping"
    } else {
        mmap = ""
    }

    """
    kraken2 \
      --threads ${task.cpus} \
      ${mmap} \
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


if ( params.merge ) {
    /*
     * Brian Bushnell recommends loose for insert size estimation, and vstrict
     * for assembly.
     */

    /*
     * Merge or stitch the read pairs.
     */
    process mergePairs {
        label "bbmap"
        label "big_task"

        publishDir "${params.outdir}/merged_reads"

        tag { ${base_name} }

        input:
        set val(base_name),
            val(fwd_name), val(rev_name),
            file("fwd"), file("rev") from joined4MergePairs

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
          loose \
          rem \
          k=62 \
          extend2=50 \
          ecct
        """
    }
}
