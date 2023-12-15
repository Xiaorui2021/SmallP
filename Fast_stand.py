"""
Objective:
1. Deduplication of the database.
2. Simplification of database titles and output of fasta files.
"""

import time
import os
import shutil

#Deduplication of the protein sequences
def pro_der(file_in,file_out):
    """
    The object is to remove duplicates from the protein Fasta sequences
    file_in: Input fasta file
    file_out: Output fasta file
    """
    from Bio import SeqIO
    seen = []
    records = []
    i = 0
    for record in SeqIO.parse(file_in, "fasta"):
        if str(record.seq) not in seen:
            seen.append(str(record.seq))
            records.append(record)
            i += 1
    # writing to a fasta file
    SeqIO.write(records, file_out, "fasta")
    file_out.close()

#Standardizing the headers of fasta sequences of proteins
def title_refine(file_in, file_out,Spe):
    """
    The object is to standardize the headers of protein Fasta sequences
    :param file_in: Input fasta file
    :param file_out: Output fasta file
    :param spe: Species
    """
    fa_Con = file_in.read()
    file_in.close()
    every_fas = fa_Con.split(">")
    for i in every_fas:
        if i != "":
            t_start = i.index(" ")
            s_start = i.index("\n")
            title_con=i[:t_start]
            seq_con = i[s_start:].replace("*","")
            file_out.write(">" +"tr|"+ title_con + "|" +" "+ "OS="+str(Spe) + " " + seq_con)
    file_out.close()

start_time = time.time()
file_path= r'C:\Users\bqxin\Desktop\SmallP\fasta data'
file_name = os.listdir(file_path) # Read all file names in the directory

#step1:Deduplication of the database
print("step1:Deduplication of the database")
file_name = [i for i in file_name if '.fasta' in i] # Keep only files with names containing ".fasta".
file_name.sort(key=lambda x: str(''.join(filter(str.isalpha, x)))) #Arrange file names in alphabetical order
for file in file_name:
    file_in = open(file_path+"/"+file, 'r')
    start = file.index(".")
    title = file[:start]
    file_out = open(file_path + "/" + title + "_rep.fasta", 'w')
    pro_der(file_in, file_out)

#step2:Compacting database titles
print("step2:Compact database titles")
file_name4 = os.listdir(file_path)
file_name4 = [i for i in file_name4 if '_rep.fasta' in i]
file_name4.sort(key=lambda x: str(''.join(filter(str.isalpha, x))))
for file in file_name4:
    file_in = open(file_path+"/"+file, 'r')
    start = file.index("_")
    title = file[:start]
    file_out = open(file_path + "/" + title + "_final.fasta", 'w')
    Spe=str(title)
    title_refine(file_in, file_out,Spe)

# step3: Delete excess files
print("Step3: Delete excess files")
fname = os.listdir(file_path)
fname = [i for i in fname if '_rep' in i]
for i in fname:
    os.remove(file_path + "/" + i)

# Computational time consumption
end_time = time.time()
escape = end_time - start_time
print('The code run {:.0f}m {:.0f}s'.format(escape // 60, escape % 60))
#————————————————
"""
Qirui Bi
Shanghai
December 7, 2023
"""
