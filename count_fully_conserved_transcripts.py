#!/usr/bin/env python

import sys

n_male_samples = 145
n_female_samples = 147


transcript_valid_for_male = set()
transcript_valid_for_female = set()
transcript_invalid = set()
transcript_chromosome = {}

for l in sys.stdin:
	parts = l.strip().split("\t")
	if parts[0] == "Chromosome": continue
	transcript = parts[2]
	female_count = int(parts[4])
	male_count = int(parts[5])
	alleleset = parts[3]
	chromosome = parts[0]
	transcript_chromosome[transcript] = chromosome
	valid = False
	if chromosome == "chrX":
		if male_count == n_male_samples and female_count == 0 and alleleset == "ref":
			transcript_valid_for_male.add(transcript)
			valid = True
		if male_count == 0 and female_count == n_female_samples and alleleset == "ref+ref":
			transcript_valid_for_female.add(transcript)
			valid = True
	elif chromosome == "chrY":
		if male_count == n_male_samples and female_count == 0 and alleleset == "ref":
			transcript_valid_for_male.add(transcript)
			valid = True
		if male_count == 0 and female_count == n_female_samples and alleleset == "missing":
			transcript_valid_for_female.add(transcript)
			valid = True
	else:
		if male_count == n_male_samples and female_count == n_female_samples and alleleset == "ref+ref":
			transcript_valid_for_male.add(transcript)
			transcript_valid_for_female.add(transcript)
			valid = True
	if not valid:
		transcript_invalid.add(transcript)

count_valids_in_chrX = 0
count_valids_in_chrY = 0
count_valids_in_autosome = 0

for x in transcript_valid_for_male:
	if x not in transcript_valid_for_female: continue
	if x in transcript_invalid: continue
	if transcript_chromosome[x] == "chrX":
		count_valids_in_chrX += 1
	elif transcript_chromosome[x] == "chrY":
		count_valids_in_chrY += 1
	else:
		count_valids_in_autosome += 1

print(f"{count_valids_in_autosome} fully conserved in autosomes")
print(f"{count_valids_in_chrX} fully conserved in chrX")
print(f"{count_valids_in_chrY} fully conserved in chrY")

