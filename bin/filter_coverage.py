#!/usr/bin/env python3

import sys
import argparse
from collections import namedtuple
from collections import defaultdict
from Bio import SeqIO
from Bio.Alphabet import generic_dna

__version__ = "0.0.1"


Coord = namedtuple(
    "Coord",
    [
        "ref_start",
        "ref_end",
        "query_start",
        "query_end",
        "ref_aln_len",
        "query_aln_len",
        "pid",
        "ref_len",
        "query_len",
        "ref_coverage",
        "query_coverage",
        "ref",
        "query",
    ]
)


Filtered = namedtuple(
    "Filtered",
    [
        "query",
        "ref",
        "query_len",
        "ref_len",
        "cov",
        "pid",
        "nmatches",
        "is_filtered",
        "reason"
    ]
)


def parse_coords(handle):

    columns = [
        ("ref_start", str),
        ("ref_end", str),
        ("query_start", int),
        ("query_end", int),
        ("ref_aln_len", int),
        ("query_aln_len", int),
        ("pid", float),
        ("ref_len", int),
        ("query_len", int),
        ("ref_coverage", float),
        ("query_coverage", float),
        ("ref", str),
        ("query", str),
    ]

    for line in handle:
        sline = line.strip().split("\t")
        assert len(sline) == len(columns)

        row = {}
        for (col, fn), val in zip(columns, sline):
            row[col] = fn(val)

        yield Coord(**row)
    return


def basic_filter_coord(coord, cov_thres):
    if coord.ref == coord.query:
        return False
    else:
        return (coord.query_coverage / 100) >= cov_thres


def grouped_filter(coords, pid_thres, max_length):
    grouped = defaultdict(list)
    out = list()

    for coord in coords:
        grouped[coord.query].append(coord)

    for group in grouped:
        grouped[group].sort(
            key=lambda x: (x.pid / 100) * (x.query_coverage / 100),
            reverse=True
        )

    for group in grouped:
        nmatches = len(grouped[group])

        first = grouped[group][0]
        too_long = first.query_len > max_length

        if too_long:
            out.append(Filtered(
                query=group,
                ref=first.ref,
                query_len=first.query_len,
                ref_len=first.ref_len,
                cov=first.query_coverage,
                pid=first.pid,
                nmatches=nmatches,
                is_filtered=False,
                reason="too long"
            ))
            continue

        best = None
        for coord in grouped[group]:
            if (coord.pid / 100) >= pid_thres:
                best = coord
                break

        if best is None:
            best = first
            is_filtered = False
            reason = "pid too low"
        else:
            is_filtered = True
            reason = ""

        out.append(Filtered(
            query=group,
            ref=best.ref,
            query_len=best.query_len,
            ref_len=best.ref_len,
            cov=best.query_coverage,
            pid=best.pid,
            nmatches=nmatches,
            is_filtered=is_filtered,
            reason=reason
        ))
    return out


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Renames fasta sequences using a sane format.",
    )

    parser.add_argument(
        "-i", "--input",
        type=argparse.FileType('r'),
        default=sys.stdin,
        help="Input coords file. Use '-' for stdin (default)."
    )

    parser.add_argument(
        "-f", "--fasta",
        type=argparse.FileType('r'),
        required=True,
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
        "-c", "--coverage",
        type=float,
        default=0.95,
        help=(
            "The coverage of the shorter sequence reqiured."
        )
    )

    parser.add_argument(
        "-p", "--pid",
        type=float,
        default=0.9,
        help=(
            "The pid of the alignment reqiured."
        )
    )

    parser.add_argument(
        "-x", "--max-length",
        type=int,
        default=10000,
        help=(
            "The maximum contig length to filter out."
        )
    )

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__
    )

    parsed = parser.parse_args(args)
    return parsed


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    seqs = SeqIO.parse(
        args.fasta,
        format="fasta",
        alphabet=generic_dna
    )
    seqs = list(seqs)
    print(seqs)

    coords = parse_coords(args.input)
    coords = (c for c in coords if basic_filter_coord(c, args.coverage))
    matches = grouped_filter(coords, args.pid, args.max_length)

    to_exclude = {c.query for c in matches if c.is_filtered}

    filtered_seqs = [s for s in seqs if s.id not in to_exclude]

    SeqIO.write(filtered_seqs, args.output, format="fasta")

    columns = [
        "query",
        "ref",
        "query_len",
        "ref_len",
        "cov",
        "pid",
        "nmatches",
        "is_filtered",
        "reason"
    ]

    print("\t".join(columns), file=args.table)
    for match in matches:
        print(
            "\t".join([str(getattr(match, col)) for col in columns]),
            file=args.table
        )
    return


if __name__ == "__main__":
    main()
