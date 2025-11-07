#!/usr/bin/env python

import sys

sample_order_file = sys.argv[1]

# exported sample allelesets from stdin
# tsv to stdout

def count_alleles(isoform_coverage):
	singletons = 0
	multiples = 0
	for transcript in isoform_coverage:
		for isoform in isoform_coverage[transcript]:
			if isoform_coverage[transcript][isoform] == 1:
				singletons += 1
			else:
				multiples += 1
	return singletons, multiples

sample_order = []
with open(sample_order_file) as f:
	for l in f:
		sample_order.append(l.strip())

sample_transcript_alleles = {}
for l in sys.stdin:
	parts = l.strip().split("\t")
	if parts[0] == "Chromosome": continue
	sample = parts[3]
	transcript = parts[2]
	alleleset = parts[4]
	if alleleset == "missing":
		alleleset = []
	else:
		alleleset = alleleset.split("+")
	if sample not in sample_transcript_alleles: sample_transcript_alleles[sample] = []
	sample_transcript_alleles[sample].append((transcript, alleleset))

print("n_samples\tadded_sample\tsingletons\tmultiples")
sample_num = 0
isoform_coverage = {}
for sample in sample_order:
	for transcript, alleleset in sample_transcript_alleles[sample]:
		if transcript not in isoform_coverage: isoform_coverage[transcript] = {}
		for isoform in alleleset:
			if isoform not in isoform_coverage[transcript]: isoform_coverage[transcript][isoform] = 0
			isoform_coverage[transcript][isoform] += 1
	singleton_alleles, multiton_alleles = count_alleles(isoform_coverage)
	sample_num += 1
	print(f"{sample_num}\t{sample}\t{singleton_alleles}\t{multiton_alleles}")
