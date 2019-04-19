#!/usr/bin/env python3

import re
import sys
import argparse
from os.path import split as psplit
from os.path import splitext
from math import ceil, log10
from Bio import SeqIO
from Bio.Alphabet import generic_dna


__version__ = "0.0.1"


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Renames fasta sequences using a sane format.",
    )

    parser.add_argument(
        "-i", "--input",
        type=argparse.FileType('r'),
        default=sys.stdin,
        help="Input fasta file. Use '-' for stdin (default)."
    )

    parser.add_argument(
        "-o", "--output",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output fasta file. Use '-' for stdout (default)."
    )

    parser.add_argument(
        "-t", "--table",
        type=argparse.FileType('w'),
        default=None,
        help=(
            "Write a tab separated file mapping the original names to the "
            "new ones."
        )
    )

    parser.add_argument(
        "-g", "--genome",
        help="A genome name, default is filename without extensions",
        default=None,
    )

    parser.add_argument(
        "-c", "--chromosome",
        default="contig",
        help=(
            "A string specifying the prefix for contigs/chromosomes, which "
            "will then have a padded number added to the end in descending "
            "sorted order of length. E.G. 'contig_' would give contig_001, "
            "contig_002 ..., assuming the total number of contigs was < 1000."
        ),
    )

    parser.add_argument(
        "-s", "--sep",
        default=".",
        help="The string to separate the genome and chromosome names by.",
    )

    parser.add_argument(
        "-p", "--spades",
        default=False,
        action="store_true",
        help=(
            "Special option for spades input, to save coverage info into "
            "description field."
        )
    )

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__
    )

    parsed = parser.parse_args(args)
    if parsed.genome is None:
        parsed.genome = psplit(splitext(parsed.input.name)[0])[1]

    return parsed


def split_spades(string):
    match = re.search(r"cov_(\d+\.?\d*)", string)
    if len(match.groups()) == 1:
        return float(match.groups()[0])
    else:
        return None


def format_description(id_, description, length, is_spades):
    template = "{}={}"
    output = []

    output.append(template.format("length", length))
    if is_spades:
        cov = split_spades(id_)
        if cov is not None:
            output.append(template.format("coverage", cov))

    if description is not None and description != "" and description != id_:
        output.append(template.format("description", description))

    return ";".join(output)


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    seqs = SeqIO.parse(args.input, format="fasta", alphabet=generic_dna)
    seqs = sorted(seqs, key=lambda s: len(s.seq), reverse=True)

    padding = ceil(log10(len(seqs) + 1))
    template = f"{args.genome}{args.sep}{args.chromosome}{{:0>{padding}}}"

    output_seqs = []
    output_tsv = []
    for i, seq in enumerate(seqs, 1):
        new_id = template.format(i)
        original_id = seq.id
        original_description = seq.description
        length = len(seq.seq)

        seq.id = new_id
        seq.name = new_id
        seq.description = format_description(
            original_id,
            original_description,
            length,
            args.spades
        )

        output_seqs.append(seq)
        output_tsv.append((args.genome, new_id,
                           original_id, original_description))

    SeqIO.write(output_seqs, args.output, format="fasta")

    if args.table is not None:
        header = ["genome", "id", "original_id", "original_description"]
        print("\t".join(header), file=args.table)

        for line in output_tsv:
            print("\t".join(line), file=args.table)

    return


if __name__ == "__main__":
    main()
