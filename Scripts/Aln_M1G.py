'''
This script performs the alignment against the mature genome allowing 1 mismatch in the seed.
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
dir_path.split('/')
path=dir_path.split('/')
path=path[:-1]
refgen_path=('/').join(path)+'/Reference_Genomes/'
sample=sys.argv[1]
sample_name=sample.split('.')[0]
func=os.path.abspath(os.getcwd())


if os.path.exists('../Results/'+sample_name):
	if not os.path.exists('../Results/'+sample_name+'/Alignment_M1G'):
			os.system('mkdir'+' '+'../Results/'+sample_name+'/Alignment_M1G')
	else:
		print ('There is a Results folder for the mature genome alignment already created for this sample ! Maybe the analysis has already been done (Cheek the Results folder)')
else:
	print ('The results folder for this sample has not been created so the first step has not been done ! ')
	raise SystemExit

os.chdir('../Results/'+sample_name+'/Alignment_M1G/')


# Alignment
print ('4-Performing the alignment against the mature genome (1 mismatch in the seed)')
command = "bowtie2 --local -p " + str(nucleos)
os.system(command+ ' -N 1 -x '+refgen_path+'families_tRNA_refgenome '+'../Alignment_PG/'+sample_name+'_unmapped_PGloc.fastq'+' --un '+sample_name+'_unmapped_MG_1M_loc.fastq'+' | samtools view -bSF4 - > '+sample_name+'_MG_1M_loc_mapped.bam')

# Processing files (sorting and indexing)
os.system('samtools sort '+sample_name+'_MG_1M_loc_mapped.bam'+ ' -o ' +sample_name+'_MG_1M_loc_mapped_sort.bam')
os.system('samtools index '+sample_name+'_MG_1M_loc_mapped_sort.bam')   

# Removing soft clipped bases
os.system('python3 '+func+'/rm_soft_clipped_bases.py '+sample_name+' '+'_MG_1M_loc_mapped_sort')

# Processing files (sorting and indexing)
os.system('samtools sort '+sample_name+'_MG_1M_loc_mapped_sort_soft_clipped_rm.bam'+ ' -o ' +sample_name+'_MG_1M_loc_mapped_sort_soft_clipped_rm_sort.bam')
os.system('samtools index '+sample_name+'_MG_1M_loc_mapped_sort_soft_clipped_rm_sort.bam') 

# Marge the files for mature tRNA
os.system('samtools merge '+'../Final_results/'+sample_name+'all_mature.bam ../Alignment_M1G/'+sample_name+'_MG_1M_loc_mapped_sort_soft_clipped_rm_sort.bam '+'../Alignment_MG/'+sample_name+'_MGloc_mapped_sort_soft_clipped_rm_sort.bam')

# Processing files (sorting and indexing)
os.system('samtools sort ../Final_results/'+sample_name+'all_mature.bam '+ ' -o ' +'../Final_results/'+sample_name+'all_mature_sort.bam')
os.system('samtools index '+'../Final_results/'+sample_name+'all_mature_sort.bam')  

