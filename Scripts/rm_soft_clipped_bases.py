#### This script removes soft clipped bases without removing the CCA tail present in some reads as soft clipped bases #####

import argparse
import subprocess
import os 
import sys
import glob
import pysam
import re



sample_name=sys.argv[1]
file=sys.argv[2]

bamfile = pysam.AlignmentFile(sample_name+file+'.bam', "rb")
out_bam = pysam.AlignmentFile(sample_name+file+'_soft_clipped_rm.bam', "wb", template=bamfile)

for read in bamfile.fetch():
	sequence=str(read.query_sequence)
	cigar_str=read.cigarstring

	if 'S' in cigar_str:

		if str(read.cigar[0])[1:-1].split(',')[0] == '4':

			soft_clipp=str(read.cigar[0])[1:-1].split(',')[1].replace(' ','')

			if soft_clipp == '3' and read.query_sequence[:3] == 'TGG':
				a=0
			else:
				soft_clipp_seq=str(sequence[:int(soft_clipp)])
				if soft_clipp_seq[-3:] == 'TGG':
					soft_clipp_cca=int(soft_clipp)-3
					cigar_soft_clipp_rm=cigar_str[(len(soft_clipp)+1):]
					cigar_soft_clipp_rm='3S'+str(cigar_soft_clipp_rm)
					read.cigarstring=cigar_soft_clipp_rm
					qual=read.query_qualities[(int(soft_clipp)-3):]
					read.query_sequence=str(sequence[(int(soft_clipp)-3):])
					read.query_qualities=qual
					cigar_str=read.cigarstring
					sequence=str(read.query_sequence)



				else:
					cigar_soft_clipp_rm=cigar_str[(len(soft_clipp)+1):]
					qual=read.query_qualities[int(soft_clipp):]
					read.cigarstring=cigar_soft_clipp_rm
					read.query_sequence=sequence[int(soft_clipp):]
					read.query_qualities=qual

					cigar_str=read.cigarstring
					sequence=str(read.query_sequence)




		if str(read.cigar[-1])[1:-1].split(',')[0] == '4':
			
			soft_clipp=str(read.cigar[-1])[1:-1].split(',')[1].replace(' ','')

			if soft_clipp == '3' and read.query_sequence[-3:] == 'CCA':

				out_bam.write(read)

			else:
				soft_clipp_seq=str(sequence[-int(soft_clipp):])
				if soft_clipp_seq[:3] == 'CCA':

					soft_clipp_cca=int(soft_clipp)-3
					cigar_soft_clipp_rm=cigar_str[:-(len(soft_clipp)+1)]
					cigar_soft_clipp_rm=str(cigar_soft_clipp_rm)+'3S'
					read.cigarstring=cigar_soft_clipp_rm
					qual=read.query_qualities[:-(int(soft_clipp)-3)]
					read.query_sequence=str(sequence[:-(int(soft_clipp)-3)])
					read.query_qualities=qual
					out_bam.write(read)

				else:
					cigar_soft_clipp_rm=cigar_str[:-(len(soft_clipp)+1)]
					read.cigarstring=cigar_soft_clipp_rm
					qual=read.query_qualities[:-int(soft_clipp)]
					read.query_sequence=str(sequence[:-int(soft_clipp)])
					read.query_qualities=qual
					out_bam.write(read)
		else:
			out_bam.write(read)
	else:
		out_bam.write(read)

out_bam.close()
bamfile.close()

