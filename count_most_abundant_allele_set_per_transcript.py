#!/usr/bin/env python

import sys

# allele sets from stdin
# counts tsv to stdout

highest_count = {}
for l in sys.stdin:
	parts = l.strip().split("\t")
	if parts[0] == "Chromosome": continue
	count = int(parts[4])
	transcript = parts[2]
	if transcript not in highest_count: highest_count[transcript] = 0
	highest_count[transcript] = max(highest_count[transcript], count)

for transcript in highest_count:
	print(transcript + "\t" + str(highest_count[transcript]))
