
'''
GUI Interface Integrated pipeline for the analisis of tRNA-Seq Datasets (IP-tRNA-dat)
'''

# Built-in/Generic Imports
import sys
import os
import platform
import subprocess

# Modules
o_sys=platform.system()

if "tkinter" not in sys.modules:
    if "Darwin" in o_sys:
        os.system("pip install tk")
    else:
        os.system("sudo apt-get install python3-tk")

import tkinter 
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
import tkinter.font as tkFont
from tkinter import messagebox as mb


# GUI Interface
app = tkinter.Tk()
app.title("tRNA Pipeline")
app.geometry('600x600+100+100')
tab_parent = ttk.Notebook(app)

tab1 = tkinter.Frame(tab_parent)
tab2 = tkinter.Frame(tab_parent)
tab2.config(bg="#F7F7F7")
tab1.config(bg="#F7F7F7")
tab_parent.add(tab1, text="IP-tRNA-dat")
tab_parent.pack(expand=1, fill='both')
lbl=tkinter.Label(app, text="")


# Functions definitions

def retrieve_SRR():
    '''
    This function recieves an SRR accession number and with string methods 
    creates the ftp url where the file is, and dowloads it. It uses different 
    methods depending on the operating system.
    ''' 
    lbl.config(text="Downloading your selected fastq.gz file.")
    SRR = text_Widget.get("1.0",'end-1c')
    first_digits=SRR[0:6]
    last_digits='00'+SRR[-1:]
    ftp_link1='ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+first_digits+'/'+last_digits+'/'+SRR+'/'+SRR+'.fastq.gz'
    ftp_link2 ="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"+first_digits+"/"+SRR+"/"+SRR+".fastq.gz"

    if 'Linux' in o_sys or 'Darwin' in o_sys:
        os.system('wget -P ' + app.fastqFolder + " " + ftp_link1)
        if SRR+".fastq.gz" not in os.listdir(app.fastqFolder):
            os.system('wget -P ' + app.fastqFolder + " " + ftp_link2)
            
        os.system("gunzip " + app.fastqFolder + "/" + SRR + ".fastq.gz")
    elif 'Windows' in o_sys:
        os.chdir(app.fastqFolder)  # This can fail!
        
        os.system('curl.exe -O '+ftp_link1)
        if SRR + ".fastq.gz" not in os.listdir(app.fastqFolder):
            os.system('curl.exe -O '+ftp_link2)
        os.system('wsl gzip -d '+SRR+'.fastq.gz')
        #os.system('del '+SRR+'.fastq.gz')
        
    
    if SRR+".fastq" in os.listdir(app.fastqFolder):
        mb.showinfo("Message", "The sample was downloaded correctly.")
    else:
        mb.showerror("Warning", "The sample was not dowloaded correctly, check the ID.")
    
    os.chdir(app.scriptsFolder)

def download_Genome():
    '''
    This function calls the scripts to download the Human genome and move 
    the files where they need to be.
    '''
    lbl.config(text="Downloading and extracting reference genome. This takes a while and freezes the app, don't close it!")
    if 'Linux' in o_sys or 'Darwin' in o_sys:
        os.system('bash download_genome_index.sh')
    elif 'Windows' in o_sys:
        os.system('wsl python3 download_genome_index.py')
    
    if ["genome.1.bt2","genome.rev.1.bt2","genome.2.bt2", "genome.rev.2.bt2",
        "genome.4.bt2"] in os.listdir(app.sourceFolder()+"/Reference_Genomes"):
        mb.showinfo("Message", "The genome was downloaded correctly.")
    else:
        mb.error("Warning", "The genome was not downloaded correctly.")
    
    
def install_programs():
    '''
    This function calls the script used to setup anaconda and the required packages.
    '''
    if 'Linux' in o_sys or 'Darwin' in o_sys:
        os.system('bash requirements.sh')
    if 'Windows' in o_sys:
        os.system("wsl python3 requirements.py")
        
    mb.showinfo("Message", "Dependencies downloaded. You can continue.")


def countsAnalysis():
    '''
    This function calls the R scripts used to create sequencing counts charts.
    '''
    os.chdir(app.scriptsFolder)
    if 'Linux' in o_sys or 'Darwin' in o_sys:
        subprocess.call(["Rscript", "--vanilla", "PLOTS_Counts_Total.r"])
        subprocess.call(["Rscript", "--vanilla", "PLOTS_Prop_Mature_Precursor.r"])
    if 'Windows' in o_sys:
        os.system('"C:\Program Files\R\R-4.0.3\bin\Rscript.exe" PLOTS_Counts_Total.r')
        os.system('"C:\Program Files\R\R-4.0.3\bin\Rscript.exe" PLOTS_Prop_Mature_Precursor.r')

def modAnalysis():
    '''
    This function calls the R scripts used to create modification ration and 
    coverage charts.
    '''
    os.chdir(app.scriptsFolder)
    if 'Linux' in o_sys or 'Darwin' in o_sys:
        subprocess.call(["Rscript", "--vanilla", "PLOTS_Modifications_Analisis.r"])
    if 'Windows' in o_sys:
        file = open("PLOTS_Modifications_Analisis.bat","w")
        file.write("Rscript.exe PLOTS_Modifications_Analisis.r %*")
        file.close()
        subprocess.call("PLOTS_Modifications_Analisis")

        
   
def run_alignment():
    '''
    This function reads the file with the samples in order to perform the tRNA 
    alignment pipeline.
    '''
    os.chdir(app.fastqFolder)
    file = open("sample_data.txt", "r")
    sample_data = file.readlines()
    sample_data = sample_data[1:]

    for sample in sample_data:
        sample = sample.split("\t")[0]
        sample = sample + ".fastq"
        alignment(sample)
    
    os.chdir(app.scriptsFolder)
    r_preparation()
    
    mb.showinfo("Message", "Alignment done!")
    

def alignment(sample):
    '''
    This function calls each step of the tRNA alignment pipeline.
    '''
    os.chdir(app.scriptsFolder)
    lbl.config(text="Starting the pipeline.")
    lbl.config(text="Checking for missing modules and installing them if missing.")
    os.system('python3 modules.py')
    lbl.config(text="Done.")
    lbl.config(text="Starting the whole genome alignment. This takes time.")
    if 'Linux' in o_sys or 'Darwin' in o_sys:
        os.system('python3 Aln_WG.py '+ sample + ' ' + app.fastqFolder)
        lbl.config(text="Done!")
        lbl.config(text="Aligning versus the mature genome.")
        os.system('python3 Aln_MG.py '+ sample)
        lbl.config(text="Done!")
        lbl.config(text="Aligning versus the precursor genome.")
        os.system('python3 Aln_PG.py '+ sample)
        lbl.config(text="Done!")
        lbl.config(text="Aligning versus the mature genome with one mismatch in the seed.")
        os.system('python3 Aln_M1G.py '+ sample)
        lbl.config(text="Done!")
        lbl.config(text="Obtaining the final counts.")
        os.system('python3 Obtain_counts.py '+ sample)
        lbl.config(text="Done!")
        lbl.config(text="Doing the pileup.")
        os.system('python3 pileup_mod.py '+ sample)
        lbl.config(text="Done!")
        
    if 'Windows' in o_sys:
        subprocess.check_call(['wsl','python3', 'Aln_WG.py ',sample,' ', app.fastqFolder])
        lbl.config(text="Done!")
        lbl.config(text="Aligning versus the mature genome.")
        subprocess.check_call(['wsl','python3', 'Aln_MG.py ',sample])
        lbl.config(text="Done!")
        lbl.config(text="Aligning versus the precursor genome.")
        subprocess.check_call(['wsl','python3', 'Aln_PG.py ',sample])
        lbl.config(text="Done!")
        lbl.config(text="Aligning versus the mature genome with one mismatch in the seed.")
        subprocess.check_call(['wsl','python3', 'Aln_M1G.py ',sample])
        lbl.config(text="Done!")
        lbl.config(text="Obtaining the final counts.")
        subprocess.check_call(['wsl','python3', 'Obtain_counts.py ',sample])
        lbl.config(text="Done!")
        lbl.config(text="Doing the pileup.")
        subprocess.check_call(['wsl','python3','pileup_mod.py ',sample])
        lbl.config(text="Done!")



def r_preparation():
    '''
    This script takes all the counts step result files and joins them
    in a single file in order to analyse them in R
    '''
    os.chdir(app.scriptsFolder)
    os.chdir("..")
    if "R_files" not in os.listdir(os.getcwd()+"/Results"):    
        os.mkdir(os.getcwd()+"/Results/R_files")
    
    os.chdir(os.getcwd()+"/Results/R_files")
    if "Base_calling" not in os.listdir(os.getcwd()):    
        os.mkdir("Base_calling")
    if "Counts" not in os.listdir(os.getcwd()):
        os.mkdir("Counts")
    if "Mapping_quality" not in os.listdir(os.getcwd()):
        os.mkdir("Mapping_quality")
        
    file = open("../../Fastq_downloaded/sample_data.txt", "r")
    samples = []
    for line in file.readlines():
        samples.append(line.split("\t")[0])
    file.close()
    file = open("../../Fastq_downloaded/sample_data.txt", "r")
    if file.readlines()[-1] != "\n":
        file.close()
        file = open("../../Fastq_downloaded/sample_data.txt", "a")    
        file.write("\n")
        file.close()
    file.close()
    samples = samples[1:]  # Remove the header.

    #sample_data.write("\n")
    #sample_data.close()

    if "Linux" in o_sys or "Darwin" in o_sys:
        os.chdir("..")
        for sample in samples:
            path_to_file1 = os.getcwd()+"/"+sample+"/Counts/"+sample+\
                "all_total.txt" 
            path_to_file2 = os.getcwd()+"/"+sample+"/Counts/"+sample+\
                "all_mature_sort.txt"                 
            path_to_file3 = os.getcwd()+"/"+sample+"/Base_calling/"+sample+ \
                "_all_base_calling_by_pos_CORRECT_OK.txt" 
            path_to_file4 = os.getcwd()+"/"+sample+"/Alignment_WG/"+sample+\
                "_WGloc_only_trna.bam"
            path_to_file5 = os.getcwd()+"/"+sample+"/Alignment_MG/"+sample+\
                "_MGloc_mapped.bam"
            path_to_file6 = os.getcwd()+"/"+sample+"/Alignment_PG/"+sample+\
                "_PGloc_mapped.bam"
            path_to_file7 = os.getcwd()+"/"+sample+"/Alignment_M1G/"+sample+\
                "_MG_1M_loc_mapped.bam"            

            cmd1 = "cp " + path_to_file1 + " " + os.getcwd() + "/R_files" + \
                "/Counts"
            cmd2 = "cp " + path_to_file2 + " " + os.getcwd() + "/R_files" + \
                "/Counts"
            cmd3 = "cp " + path_to_file3 + " " + os.getcwd() + "/R_files" + \
                "/Base_calling"
            cmd4 = "cp " + path_to_file4 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
            cmd5 = "cp " + path_to_file5 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
            cmd6 = "cp " + path_to_file6 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
            cmd7 = "cp " + path_to_file7 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
                
            os.system(cmd1)
            os.system(cmd2)
            os.system(cmd3)
            os.system(cmd4)
            os.system(cmd5)
            os.system(cmd6)
            os.system(cmd7)
        os.chdir(app.scriptsFolder)
        subprocess.call(["Rscript", "--vanilla", "Join_results.r"])
        
    if "Windows" in o_sys:  # This is gonna fail.
        os.chdir("..")
        for sample in samples:
            path_to_file1 = os.getcwd()+"\ "+sample+"\Counts\ "+sample+\
                "all_total.txt" 
            path_to_file2 = os.getcwd()+"\ "+sample+"\Counts\ "+sample+\
                "all_mature_sort.txt"    
            path_to_file3 = os.getcwd()+"\ "+sample+"\Base_calling\ "+sample+ \
                "_all_base_calling_by_pos_CORRECT_OK.txt"
            path_to_file4 = os.getcwd()+"\ "+sample+"\Alignment_WG\ "+sample+\
                "_WGloc_only_trna.bam"
            path_to_file5 = os.getcwd()+"\ "+sample+"\Alignment_MG\ "+sample+\
                "_MGloc_mapped.bam"
            path_to_file6 = os.getcwd()+"\ "+sample+"\Alignment_PG\ "+sample+\
                "_PGloc_mapped.bam"
            path_to_file7 = os.getcwd()+"\ "+sample+"/Alignment_M1G\ "+sample+\
                "_MG_1M_loc_mapped.bam"  
                
            cmd1 = "copy " + path_to_file1 + " " + os.getcwd() + "\R_files" + \
                "\Counts"
            cmd2 = "copy " + path_to_file1 + " " + os.getcwd() + "\R_files" + \
                "\Counts"
            cmd3 = "copy " + path_to_file3 + " " + os.getcwd() + "\R_files" + \
                "\Base_calling"
            cmd4 = "copy " + path_to_file4 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
            cmd5 = "copy" + path_to_file5 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
            cmd6 = "copy " + path_to_file6 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
            cmd7 = "copy " + path_to_file7 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
                
            os.system(cmd1)
            os.system(cmd2)    
            os.system(cmd3)
            os.system(cmd4)
            os.system(cmd5)
            os.system(cmd6)
            os.system(cmd7)
        
        os.chdir(app.scriptsFolder)
        file = open("Join_results.bat","w")
        file.write("Rscript.exe Join_results.r %*")
        file.close()
        subprocess.call("Join_results")
    
    

def remove_alignments():
    """
    This function remove the aligments files.
    """
    file = open("../Fastq_downloaded/sample_data.txt", "r")
    samples = []
    for line in file.readlines():
        samples.append(line.split("\t")[0])
    samples = samples[1:]  # Remove the header.
    
    os.chdir("../Results")
    for sample in samples: 
        if sample in os.listdir():
            if "Linux" in o_sys or "Darwin" in o_sys:
                os.system("rm -r " + sample)
            
            if "Windows" in o_sys:
                os.system("rd /s /q" + sample)
                
        else:
            pass
        
    os.chdir(app.scriptsFolder)
    

def r_analysis():
    '''
    This function call the count and modification analysis Rscripts. And launches an Rscript 
    for creating heatmaps depending on the modification ratio and the final report.
    '''
    file = open("../Fastq_downloaded/sample_data.txt", "r")
    if file.readlines()[-1] != "\n":
        file.close()
        file = open("../Fastq_downloaded/sample_data.txt", "a")    
        file.write("\n")
        file.close()
    file.close()

    countsAnalysis()
    
    modAnalysis()
    
    if "Linux" in o_sys or "Darwin" in o_sys:
        subprocess.call(["Rscript", "--vanilla", "General_plots.r"])
        subprocess.call(["Rscript", "--vanilla", "Heatmaps.r"])
        subprocess.call(["Rscript", "--vanilla", "Modification_comparison_test.r"])
        subprocess.call(["Rscript", "--vanilla", "DEG_analysis.r"])
        
    if "Windows" in o_sys: 
        file = open("General_plots.bat","w")
        file.write("Rscript.exe General_plots.r %*")
        file.close()
        file2 = open("Heatmaps.bat", "w")
        file2.write("Rscript.exe Heatmaps.r %*")
        file2.close()
        subprocess.call("General_plots")
        subprocess.call("Heatmaps")
    
    mb.showinfo("Message", "Analysis done! You can check the results in the results folder.")
        

def create_text():
    '''
    This function will open a txt with the sample information, 
    user will add sample specifications
    '''
    sourceFolder=app.sourceFolder
    os.chdir(sourceFolder+"/Fastq_downloaded")
    files = os.listdir(os.getcwd())
    if "sample_data.txt" in files:
        if "Linux" in o_sys or "Darwin" in o_sys:
            os.system("rm sample_data.txt")
        if "Windows" in o_sys:
            os.system("del sample_data.txt")
        files = os.listdir(os.getcwd())
        files = [x for x in files if not x.startswith("._")]

        
    sample_data = open("sample_data.txt", "w")
    sample_data.write("ID" + "\t" + "Condition" + "\n")
    for i in range(len(files)-1):
        sample_data.write(files[i][:-6] + "\t" + "\n")
    sample_data.write(files[-1][:-6] + "\t")
    sample_data.close()
    if "Linux" in o_sys:
        os.system('gedit sample_data.txt')
    if "Darwin" in o_sys:
        os.system('open -a TextEdit sample_data.txt')
    if "Windows" in o_sys:
        os.system('notepad sample_data.txt')

    scriptsFolder=app.scriptsFolder       
    os.chdir(scriptsFolder)
        
    

# WIDGETS FOR THE GUI

app.scriptsFolder = os.getcwd()
app.sourceFolder = app.scriptsFolder[:-7]

os.chdir(app.sourceFolder)

try:
    os.mkdir("Fastq_downloaded")
except:
    pass

app.fastqFolder = app.sourceFolder + "Fastq_downloaded"
os.chdir(app.scriptsFolder)

myFont = tkFont.Font(family="Calibri", size=15)
myFont2 = tkFont.Font(family="Calibri", size=14)
myFont3 = tkFont.Font(family="Calibri", size=12)

text1=tkinter.Label(text="Set Up", width=10,height=1, bg="#F7F7F7")
text1.place(x=70, y=55)
text1.config(font='Helvetica 14 bold')

programs_Button=tkinter.Button(tab1, text='Download dependencies', width=19 , height=2,command=install_programs)
programs_Button.place(x=60, y=50)
programs_Button.config(font=myFont2)

genome_Button=tkinter.Button(tab1, text='Download Genome', width=16, height=2, command=download_Genome)
genome_Button.place(x=335, y=50)
genome_Button.config(font=myFont2)

text1=tkinter.Label(text="Select Samples", width=15, height=1,bg="#F7F7F7")
text1.place(x=75, y=165)
text1.config(font='Helvetica 14 bold')

text_Widget = tkinter.Text(height=1, width=10)

text_Widget.config(font=myFont3)
text_Widget.insert(tkinter.END, "e.g SRR44089")
text_Widget.pack()
def default(event):
    current = text_Widget.get("1.0", tkinter.END)
    if current == "e.g SRR44089\n":
        text_Widget.delete("1.0", tkinter.END)
    elif current == "\n":
        text_Widget.insert("1.0", "e.g SRR44089")

text_Widget.bind("<FocusIn>", default)
text_Widget.bind("<FocusOut>", default)
text_Widget.place(x=220, y=195)

download_Button=tkinter.Button(tab1, text='Download sample:', width=13, height=1, command=retrieve_SRR)
download_Button.place(x=55, y=160)
download_Button.config(font=myFont2)

open_text = tkinter.Button(tab1, text="Indicate sample information", width=21, height=2,
                           command=create_text)
open_text.place(x=310, y = 150)
open_text.config(font=myFont2)

text1=tkinter.Label(text="tRNA Analysis", width=15, height=1,bg="#F7F7F7")
text1.place(x=240, y=280)
text1.config(font='Helvetica 14 bold')

submit_button = tkinter.Button(tab1, text="Run alignment", width = 15, height = 2, command=run_alignment)
submit_button.place(x=200, y=280)
submit_button.config(font=myFont)

b_countsAnalysis=tkinter.Button(tab1, text='Results Report', width=15, height=2, command=r_analysis)
b_countsAnalysis.place(x=200, y= 350)
b_countsAnalysis.width=250
b_countsAnalysis.config(font=myFont)

text1=tkinter.Label(text="Gene Translation Laboratory | IRB Barcelona", width=40, height=1,bg="#F7F7F7")
text1.place(x=175, y=530)
text1.config(font='Helvetica 9 italic')

text1=tkinter.Label(text="Integrated pipeline for the analisis of tRNA-Seq Datasets", width=50, height=1,bg="#F7F7F7")
text1.place(x=145, y=550)
text1.config(font='Helvetica 9 italic')

app.mainloop()
