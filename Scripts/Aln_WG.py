'''
This script performs the alignment against the precursor genome.
'''
# Built-in/Generic Imports and Modules
import argparse
import subprocess
import os 
import sys
import glob
import re
import multiprocessing
import pysam

nucleos = multiprocessing.cpu_count()

# Sample name
sample=sys.argv[1]
folder=sys.argv[2]

print ('ANALYZING SAMPLE:',sample)

# Samples folder.
samples_path=folder

#Folder with the reference genome:
refgen_path='../Reference_Genomes/'
func=os.path.abspath(os.getcwd())



#Name of all the samples in the folder:
samples=os.listdir(samples_path)

def SAM_iterator(sam_filename):
    '''Parse function for Sam file'''
    samfile = pysam.AlignmentFile(sam_filename, "rb")
    for line in samfile:
        list_read=str(line)
        read=(list_read.split('\t'))
        yield (read)


if sample in samples:
    sample_1=sample.split('.')
    sample_name=sample_1[0]

    if sample_1[1]=='fa' or sample_1[1]=='fastq':
        #####ALIGNMENT WITH THE WHOLE GENOME#######
	
        if not os.path.exists('../Results/'+sample_name+'/Alignment_WG/'):
            os.makedirs('../Results/'+sample_name+'/Alignment_WG/')
        
        print ('1-Performing the Alignment against the whole genome:')
        command = "bowtie2 --local -p "+ str(nucleos)
        
        os.system(command + ' -N 0 -x '+refgen_path+'genome'+' '+samples_path + "/"+sample+' --un ../Results/'+sample_name+'/Alignment_WG/'+sample_name+'_unmapped_WGloc.fastq'+' | samtools view -bSF4 - > '+'../Results/'+sample_name+'/Alignment_WG/'+sample_name+'_WGloc_mapped.bam')
        
        #Set working directory.
        os.chdir('../Results/'+sample_name+'/Alignment_WG/')

        #Processing files (sorting and indexing).
        os.system('samtools sort '+sample_name+'_WGloc_mapped.bam'+ ' -o ' +sample_name+'_WGloc_mapped_sort.bam')
        os.system('samtools index '+sample_name+'_WGloc_mapped_sort.bam')   
		


        print ('1.2-Filtering reads. Obtain reads that are mapping to tRNAs')
		
        #Obtain the reads that correspond to tRNA only.
        os.system('samtools view -b -L '+'../../../Reference_Genomes/info/tRNA_only.bed'+' '+sample_name+'_WGloc_mapped.bam > '+sample_name+'_WGloc_only_trna.bam')
		
        #Processing files (sorting and indexing).
        os.system('samtools sort '+sample_name+'_WGloc_only_trna.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_sort.bam')
        os.system('samtools index '+sample_name+'_WGloc_only_trna_sort.bam') 

        
        #Removing soft clipped bases
        os.system('python3 '+func+'/rm_soft_clipped_bases.py '+sample_name+' '+'_WGloc_only_trna_sort')

        #Processing files (sorting and indexing).
        os.system('samtools sort '+sample_name+'_WGloc_only_trna_sort_soft_clipped_rm.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam')
        os.system('samtools index '+sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam') 



        print ('1.3-Filtering and calssification of the reads. Classifing reads between precursor and mature sequencies')
		
        #Extract the reads that correspond to tRNA precursor only. lead_trail.bed contains the trailing regions, if a read falls in to these regions is considered a precursor. 
        os.system('samtools view -b -L '+'../../../Reference_Genomes/info/lead_trail_ok.bed '+sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam > '+sample_name+'_WGloc_only_trna_precursor.bam')
		
		#Sort the bam file.
        os.system('samtools sort '+sample_name+'_WGloc_only_trna_precursor.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_precursor_sort.bam')
		
		#Index the bam file.
        os.system('samtools index '+sample_name+'_WGloc_only_trna_precursor_sort.bam')   
		
		#Obtain the reads in sam format.
        os.system('samtools view -h -o '+sample_name+'_WGloc_only_trna_precursor_sort.sam '+sample_name+'_WGloc_only_trna_precursor_sort.bam')
		
		#Filter the precursor file to detect reads that have a deletion so they don't have the intronic region and they should be mature reads. 	
        out_id= open(sample_name+'intron_delet_mature.txt',"w")
        count_mapq=0
        all_reads=0
        fd = open(sample_name+'_WGloc_only_trna_precursor_sort.sam',"r")
        for aln in SAM_iterator(fd):
            qname=str(aln[0])
            cigar=str(aln[5])
            seq=str(aln[9])
            if 'D' in cigar:
                for num1, i_or_d, num2, m in re.findall('(\d+)([ID])(\d+)?([A-Za-z])?', cigar):
                    if int(num1) > 3:
                        out_id.write(qname+'\n')
        out_id.close()
	
	with open(sample_name+'intron_delet_mature.txt') as file2:
            first = file2.read(1)
        
	if not first:
            os.rename(sample_name+'_WGloc_only_trna_precursor_sort.bam',sample_name+'_WGloc_only_trna_with_intron_precursor_sort.bam') 
            os.system('samtools view -F 4 '+sample_name+'_WGloc_only_trna_with_intron_precursor_sort.bam '+'| cut -f1 | sort -u > '+sample_name+'_WGloc_only_trna_precursor_filtered.txt')

        else:
            #change comand line !!!! https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
            os.system('picard FilterSamReads '+'I='+sample_name+'_WGloc_only_trna_precursor_sort.bam '+'O='+sample_name+'_WGloc_only_trna_with_intron_precursor.bam READ_LIST_FILE='+sample_name+'intron_delet_mature.txt FILTER=excludeReadList >/dev/null 2>&1')
		
        	#Sort the bam file.
            os.system('samtools sort '+sample_name+'_WGloc_only_trna_with_intron_precursor.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_with_intron_precursor_sort.bam')

        	#Index the bam file.
            os.system('samtools index '+sample_name+'_WGloc_only_trna_with_intron_precursor_sort.bam') 
        	
        	#Extract the reads id of the precursor 
            os.system('samtools view -F 4 '+sample_name+'_WGloc_only_trna_with_intron_precursor_sort.bam '+'| cut -f1 | sort -u > '+sample_name+'_WGloc_only_trna_precursor_filtered.txt')

		#Extact the mature reads:
        os.system('picard FilterSamReads '+'I='+sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam '+'O='+sample_name+'_WGloc_only_trna_mature.bam READ_LIST_FILE='+sample_name+'_WGloc_only_trna_precursor_filtered.txt FILTER=excludeReadList >/dev/null 2>&1')
		
		#Sort the bam file.
        os.system('samtools sort '+sample_name+'_WGloc_only_trna_mature.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_mature_sort.bam')
		
		#Index the bam file of tRNA precursor data only:
        os.system('samtools index '+sample_name+'_WGloc_only_trna_mature_sort.bam')
		
		#Marge the mature reads +  unmmaped reads, in order to peform the aligment with the mature genome:
        os.system('bedtools bamtofastq -i '+sample_name+'_WGloc_only_trna_mature.bam -fq '+sample_name+'_WGloc_only_trna_mature.fastq')

		#Marge fastq 
        os.system('cat '+sample_name+'_WGloc_only_trna_mature.fastq '+sample_name+'_unmapped_WGloc.fastq > '+sample_name+'_WGloc_only_trna_mature_and_unmapped.fastq')

