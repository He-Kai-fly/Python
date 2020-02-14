#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#导入模块，初始传递命令
import argparse

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 fq_fa.py -i [fastq] -o [fasta]\n用于提取 fastq 中的序列并转为 fasta')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-i', '--input', metavar = '[fastq]', help = '输入文件，fastq 文件', required = True)
required.add_argument('-o', '--output', metavar = '[fasta]', help = '输出文件，fasta 文件', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

#读取 fastq 文件并提取序列转换为 fasta 文件
fasta = open(args.output, 'w')
with open(args.input, 'r') as fastq:
	for line in fastq:
		if line[0] == '@':
			line = line.strip('[@\n]')
			print('>' + line, file = fasta)
			print(fastq.readline().strip(), file = fasta)

fastq.close()
fasta.close()
