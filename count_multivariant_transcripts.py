#!/usr/bin/env python

import sys

# exported isoforms from stdin
# summary to stdout

def get_ref_isoform(transcripts_list):
	for t in transcripts_list:
		if t[0] == "ref":
			return t
	highest = transcripts_list[0]
	for i in range(1, len(transcripts_list)):
		if transcripts_list[i][1] > highest[1]:
			highest = transcripts_list[i]
	return highest

def is_multivariant(ref, alt):
	start_clip = 0
	end_clip = 0
	while start_clip < len(ref[2]) and start_clip < len(alt[2]) and ref[2][start_clip] == alt[2][start_clip]:
		start_clip += 1
	while end_clip < len(ref[2]) and end_clip < len(alt[2]) and ref[2][-1-end_clip] == alt[2][-1-end_clip]:
		end_clip += 1 
	if start_clip + end_clip >= len(ref[2]):
		return False # single indel
	if start_clip + end_clip >= len(alt[2]):
		return False # single indel
	if start_clip + end_clip + 1 == len(ref[2]) and len(ref[2]) == len(alt[2]):
		return False # single SNP
	return True # not a single indel or SNP

isoforms_per_transcript = {}
count_isoforms = 0
count_total_copies = 0
for l in sys.stdin:
	parts = l.strip().split("\t")
	if parts[0] == "Chromosome": continue
	chromosome = parts[0]
	gene = parts[1]
	transcript = parts[2]
	isoform_name = parts[3]
	coverage = int(parts[4])
	sequence = parts[5]
	if transcript not in isoforms_per_transcript: isoforms_per_transcript[transcript] = []
	isoforms_per_transcript[transcript].append((isoform_name, coverage, sequence, gene, chromosome))
	count_isoforms += 1
	count_total_copies += coverage

count_transcripts_with_multivariants = 0
count_isoforms_with_multivariants = 0
count_coverage_of_isoforms_with_multivariants = 0
count_more_than_half_has_multivariant = 0
for transcript in isoforms_per_transcript:
	ref_isoform = get_ref_isoform(isoforms_per_transcript[transcript])
	transcript_has_multivariant = False
	count_multivariants = 0
	count_not_multivariants = 0
	for isoform in isoforms_per_transcript[transcript]:
		if isoform == ref_isoform:
			count_not_multivariants += isoform[1]
			continue
		if is_multivariant(ref_isoform, isoform):
			transcript_has_multivariant = True
			count_isoforms_with_multivariants += 1
			count_coverage_of_isoforms_with_multivariants += isoform[1]
			count_multivariants += isoform[1]
		else:
			count_not_multivariants += isoform[1]
	if transcript_has_multivariant:
		count_transcripts_with_multivariants += 1
	if count_multivariants > count_not_multivariants:
		count_more_than_half_has_multivariant += 1
		print(f"majority multivariant transcript {ref_isoform[3]} {ref_isoform[4]} {transcript} {float(count_multivariants) / float(count_multivariants + count_not_multivariants)}")

print(f"{count_transcripts_with_multivariants} transcripts with multivariant")
print(f"{count_isoforms_with_multivariants} isoforms with multivariant ({float(count_isoforms_with_multivariants)/float(count_isoforms)})")
print(f"{count_coverage_of_isoforms_with_multivariants} copies with multivariant ({float(count_coverage_of_isoforms_with_multivariants)/float(count_total_copies)})")
print(f"{count_more_than_half_has_multivariant} transcripts have >50% multivariant")
