#!/usr/bin/env python

import sys

def get_copy_count(alleleset):
	if alleleset == "missing": return 0
	return alleleset.count('+') + 1

n_male_samples = 145
n_female_samples = 147

transcript_nonref_copycount = {}
transcript_chromosome = {}

for l in sys.stdin:
	parts = l.strip().split("\t")
	if parts[0] == "Chromosome": continue
	transcript = parts[2]
	if transcript not in transcript_nonref_copycount: transcript_nonref_copycount[transcript] = 0
	female_count = int(parts[4])
	male_count = int(parts[5])
	alleleset = parts[3]
	chromosome = parts[0]
	transcript_chromosome[transcript] = chromosome
	copycount = get_copy_count(alleleset)
	valid = False
	if chromosome == "chrX":
		if get_copy_count(alleleset) != 1:
			transcript_nonref_copycount[transcript] += male_count
		if get_copy_count(alleleset) != 2:
			transcript_nonref_copycount[transcript] += female_count
	elif chromosome == "chrY":
		if get_copy_count(alleleset) != 1:
			transcript_nonref_copycount[transcript] += male_count
		if get_copy_count(alleleset) != 0:
			transcript_nonref_copycount[transcript] += female_count
	else:
		if get_copy_count(alleleset) != 2:
			transcript_nonref_copycount[transcript] += male_count
			transcript_nonref_copycount[transcript] += female_count

count_variable_in_chrX = 0
count_variable_in_chrY = 0
count_variable_in_autosome = 0
for t in transcript_nonref_copycount:
	if transcript_nonref_copycount[t] <= (n_male_samples + n_female_samples) / 2:
		continue
	if transcript_chromosome[t] == "chrX":
		count_variable_in_chrX += 1
	elif transcript_chromosome[t] == "chrY":
		count_variable_in_chrY += 1
	else:
		count_variable_in_autosome += 1

print(f"{count_variable_in_autosome} majority copy count variable in autosomes")
print(f"{count_variable_in_chrX} majority copy count variable in chrX")
print(f"{count_variable_in_chrY} majority copy count variable in chrY")
