#!/usr/bin/env python

import sys

remove_contigs = set(sys.argv[1:])

printing = True
for l in sys.stdin:
	if l[0] == ">":
		contig_name = l[1:].strip()
		printing = True
		if contig_name in remove_contigs: printing = False
	if printing:
		print(l.strip())
