'''
This script performs the alignment against the mature genome.
'''

# Built-in/Generic Imports 
import os 
import sys
import subprocess
import multiprocessing

nucleos = multiprocessing.cpu_count()

# Arguments
fullCmdArguments= sys.argv
argumentList= fullCmdArguments[1:]
dir_path = os.path.dirname(os.path.realpath(__file__))
listdir=os.listdir(dir_path)
fastq=sys.argv[1]
sample_name=(fastq.split('.'))[0]
dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path.split('/')
path=dir_path.split('/')
path=path[:-1]
refgen_path=('/').join(path)+'/Reference_Genomes/'

if os.path.exists('../Results/'+sample_name):
	if not os.path.exists('../Results/'+sample_name+'/Alignment_MG'):
			os.makedirs('../Results/'+sample_name+'/Alignment_MG')
	else:
		print ('There is a Results folder for the mature genome alignment already created for this sample ! Maybe the analysis has already been done (Cheek the Results folder)')
else:
	print ('The results folder for this sample has not been created so the first step has not been done ! ')
	raise SystemExit

func=os.path.abspath(os.getcwd())

os.chdir('../Results/'+sample_name+'/Alignment_MG/')

# Alignment
print ('2-Performing the alignment against the mature genome')
command = "bowtie2 --local -p " + str(nucleos)
os.system(command + ' -N 0 -x '+refgen_path+'families_tRNA_refgenome ../Alignment_WG/'+sample_name+'_WGloc_only_trna_mature_and_unmapped.fastq'+' --un '+sample_name+'_unmapped_MGloc.fastq'+' | samtools view -bSF4 - > '+sample_name+'_MGloc_mapped.bam')

# Processing files (sorting and indexing).
os.system('samtools sort '+sample_name+'_MGloc_mapped.bam'+ ' -o ' +sample_name+'_MGloc_mapped_sort.bam')
os.system('samtools index '+sample_name+'_MGloc_mapped_sort.bam') 
 
# Removing soft clipped bases
os.system('python3 '+func+'/rm_soft_clipped_bases.py '+sample_name+' '+'_MGloc_mapped_sort')

# Processing files (sorting and indexing).
os.system('samtools sort '+sample_name+'_MGloc_mapped_sort_soft_clipped_rm.bam'+ ' -o ' +sample_name+'_MGloc_mapped_sort_soft_clipped_rm_sort.bam')
os.system('samtools index '+sample_name+'_MGloc_mapped_sort_soft_clipped_rm_sort.bam') 

	
#Join the unmapped reads with the set of precursor reads.
os.system('bedtools bamtofastq -i ../Alignment_WG/'+sample_name+'_WGloc_only_trna_with_intron_precursor_sort.bam -fq '+sample_name+'_WGloc_only_trna_with_intron_precursor_sort.fastq')
os.system('cat '+sample_name+'_WGloc_only_trna_with_intron_precursor_sort.fastq '+sample_name+'_unmapped_MGloc.fastq > '+sample_name+'_WGloc_only_trna_precursor_and_MGunmapped.fastq')

