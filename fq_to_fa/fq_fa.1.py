#!/usr/bin/env python3
# -*- coding: utf-8 -*-

fasta = open('test.fasta', 'w')
with open('test.fastq', 'r') as fastq:
	for line in fastq:
		if line[0] == '@':
			line = line.strip('[@\n]')
			print('>' + line, file = fasta)
			print(fastq.readline().strip(), file = fasta)

fastq.close()
fasta.close()
