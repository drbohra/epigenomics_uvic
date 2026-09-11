#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    description="Find the gene whose start coordinate is closest to a regulatory element."
)

parser.add_argument(
    "--input",
    required=True,
    help="Tab-separated file containing gene name and gene start coordinate"
)

parser.add_argument(
    "--start",
    required=True,
    type=int,
    help="Start coordinate of the regulatory element"
)

args = parser.parse_args()

closest_gene = None
closest_start = None
closest_distance = None

with open(args.input) as infile:
    for line in infile:
        line = line.strip()

        if not line:
            continue

        gene, start = line.split("\t")[:2]
        start = int(start)

        distance = abs(start - args.start)

        if closest_distance is None or distance < closest_distance:
            closest_gene = gene
            closest_start = start
            closest_distance = distance

print("{}\t{}\t{}".format(closest_gene, closest_start, closest_distance))
