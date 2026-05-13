#!/usr/bin/env python3
"""
Count sequences in a selection of Fasta/q files
"""
import argparse
from functools import partial
import gzip
import pathlib

GZ_EXTS = {".gz", ".gzip", ".bgz", ".bgzip"}

def open_seq(path):
    """
    Open a potentially Gzipped sequence file. Returns an open file handle in
    read-only text mode and "fasta"/"fastq" depending on the detected mode
    """
    exts = pathlib.Path(path).suffixes
    if exts[-1] in GZ_EXTS:
        _open = partial(gzip.open, mode="rt")
        exts = exts[:-1]
    else:
        _open = partial(open, mode="r")

    with _open(path) as handle:
        l = handle.readline()
        if l[0] == "@":
            mode = "fastq"
        elif l[0] == ">":
            mode = "fasta"
        else:
            raise ValueError(
                f"{path} doesn't start with an @ or >, unable to determine Fasta/Fastq status"
            )

    return _open(path), mode

def main():
    """
    Main
    """
    args = parse_args()

    pairing = ["R1", "R2"] if args.paired else ["R1"]

    if len(args.reads) != len(pairing):
        raise ValueError(f"Expected {len(pairing)} file(s), got {len(args.reads)}")

    for paired, path in zip(pairing, args.reads):
        seq_file, mode = open_seq(path)

        try:
            if mode == "fasta":
                count = sum(1 for line in seq_file if line.startswith(">"))
            elif mode == "fastq":
                count = sum(1 for _ in seq_file)

                if not count % 4 == 0:
                    raise ValueError(f"FASTQ line count is not divisible by 4: {path}")

                count = count // 4
            else:
                # Should never trigger as mode assigned in open_seq as fasta/fastq
                raise ValueError(f"None fastq/a mode for {path}")
        finally:
            seq_file.close()

        # sample stage label read file format records
        print(
            args.sample, args.stage, args.label, paired,
            pathlib.Path(path).name, mode, count, sep="\t"
        )

def arg_parser():
    """
    Construct argument parser
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sample", default="")

    parser.add_argument("--stage", default=0)

    parser.add_argument("--label", default="raw")

    parser.add_argument("--paired", action="store_true")

    parser.add_argument("reads", nargs="+")

    return parser

def parse_args():
    """
    Parse and validate script arguments
    """
    args = arg_parser().parse_args()

    return args

if __name__ == "__main__":
    main()