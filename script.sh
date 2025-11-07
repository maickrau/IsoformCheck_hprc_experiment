#!/usr/bin/sh

IsoformCheck/src/IsoformCheck stats -db db > stats.txt

IsoformCheck/src/IsoformCheck listsamples -db db > listed_samples.tsv

IsoformCheck/src/IsoformCheck exportallelesets -db db --include-gene-info > exported_allelesets.tsv

IsoformCheck/src/IsoformCheck contingencytable -db db --table sample_male_female.tsv --include-gene-info > male_female_contingency_table.tsv

IsoformCheck/src/IsoformCheck exportisoforms -db db --include-gene-info > exported_isoforms.tsv

awk 'NR==1{print "Chromosome\tGene\tTranscript\tAlleleset\tCount";}NR>1&&$4!="reference"{counts[$1 "\t" $2 "\t" $3][$5] += 1;}END{for (c in counts) { for (a in counts[c]) { print c "\t" a "\t" counts[c][a]; }}}' < exported_allelesets.tsv > alleleset_counts.tsv

awk '{counts[$1 "\t" $2 "\t" $3] -= ($5/296)*log($5/296);}END{for (key in counts) { print counts[key] "\t" key; }}' < alleleset_counts.tsv | sort -gr > most_alleleset_variable.txt

awk '$4=="missing"{print;}$4!="missing"{print $1 "\t" $2 "\t" $3 "\t" $4 "+" "\t" $5;}' < alleleset_counts.tsv | sed 's/[refABCDEFGHIJKLMNOPQRSTUVWXYZ]\++/+/g' | awk '$4!="missing"{print $1 "\t" $2 "\t" $3 "\t" length($4) "\t" $5;}$4=="missing"{print $1 "\t" $2 "\t" $3 "\t" "0" "\t" $5;}' | awk '{sums[$1 "\t" $2 "\t" $3 "\t" $4] += $5;}END{for (c in sums) {print c "\t" sums[c];}}' | awk '{entropy[$1 "\t" $2 "\t" $3] -= $5/296*log($5/296);}END{for (c in entropy) {print entropy[c] "\t" c;}}' | sort -gr > most_copy_count_variable.txt

python ./count_fully_conserved_transcripts.py < male_female_contingency_table.tsv > stats_fully_conserved.txt

python ./count_copycount_variable_transcripts.py < male_female_contingency_table.tsv > stats_copycount_variability.txt

python ./count_multivariant_transcripts.py < exported_isoforms.tsv > multivariant_stats.txt

python ./count_most_abundant_allele_set_per_transcript.py < alleleset_counts.tsv > most_abundant_per_transcript.tsv

IsoformCheck/src/IsoformCheck chisquare -db db --table merged_sample_populations.tsv --include-gene-info > chisquare_merged.tsv
IsoformCheck/src/IsoformCheck chisquare -db db --table hprc_sample_populations.tsv --include-gene-info > chisquare_hprc.tsv
IsoformCheck/src/IsoformCheck chisquare -db db --table hgsvc_sample_populations.tsv --include-gene-info > chisquare_hgsvc.tsv

cut -f 1 < listed_samples.tsv | grep -v Sample | grep -v reference | awk 'BEGIN{print "reference";} {print;}' > sample_list.txt
python ./measure_allelic_growth.py sample_list.txt < exported_allelesets.tsv > allelic_growth.tsv

Rscript plot_allelic_growth.Rscript
