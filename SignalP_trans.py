"""
Purpose: Statistical analysis of signal peptide-containing protein sequences on the basis of signal peptide prediction
Experimental content:
1. Import and read the result files of SignalP 6.0, including prediction results and region output.
2. Combine the two files according to the protein sequence titles;
3. Select the signal peptide sequences with the highest signal peptide prediction scores, and delete the sequences with lower scores
4. Output the protein sequences containing signal peptide and the mature protein fasta sequences according to the sequence scoring by SignalP 6.0.
5. Perform Cframe statistics for mature protein fasta sequences.
"""

from numpy import *
import pandas as pd
import re
import os
import time

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

start_time = time.time()
file_path1 = r'C:\Users\bqxin\Desktop\SmallP\raw data'
file_path2 = r'C:\Users\bqxin\Desktop\SmallP\fasta data'
file_path3 = r'C:\Users\bqxin\Desktop\SmallP\SignalP'
seq_screen="C"
specie="Pheretima aspergillum"

#1 Normalizing "region_output.gff3" file
print("1 Normalizing the region_output.gff3 file")
df = pd.read_table(file_path1 + "/" +  "region_output.gff3", skiprows=[0])
df.columns = ['# ID', 'Software', 'Region', 'Region_L', 'Region_R', '6', '7', '8', '9']
df=df.drop(columns=['6', '7', '8', '9'])
df.to_excel(file_path3 + "/" + "region_output_temp.xlsx", index=False)

#2 Reading the region_output.gff3 file
print("2 Reading the region_output.gff3 file")
df =pd.read_excel(file_path3 + "/" + "region_output_temp.xlsx")
col_name = df.columns.tolist()
add_list=["N_pos1", "N_pos2", "H_pos1", "H_pos2", "C_pos1", "C_pos2", "lipid_C","Seq_N","Seq_H","Seq_C","Seq_lipd", "SeqID", "SeqID_rev","rev_M","pos_M","pos_SP","sp_score","Sq_Len","C_num","C_Pos","C_frame"]
for i in range(len(add_list)):
    index = len(col_name) +i+1
    col_name.insert(index, add_list[i])
    df=df.reindex(columns=col_name)
for i in range(0,len(df["# ID"])-2):
    rate=100*i/len(df["# ID"])
    start1 = df["# ID"][i].index("-")
    start2 = df["# ID"][i].index("_P")
    SeqID = df["# ID"][i][:start1]
    SeqID_rev=df["# ID"][i][:start2]
    rev = df["# ID"][i][start1+1:start2]
    a = [m.start() for m in re.finditer("_", df["# ID"][i])]
    pos_=a[len(a)-1]+1
    mod = df["# ID"][i][pos_:]
    if df["Region"][i] == "n-region":
        df["N_pos1"][i] = df["Region_L"][i]
        df["N_pos2"][i] = df["Region_R"][i]
        df["H_pos1"][i] = df["Region_L"][i + 1]
        df["H_pos2"][i] = df["Region_R"][i + 1]
        df["C_pos1"][i] = df["Region_L"][i + 2]
        df["C_pos2"][i] = df["Region_R"][i + 2]
        df["SeqID"][i] = SeqID
        df["SeqID_rev"][i] = SeqID_rev
        df["rev_M"][i] = rev
        df["pos_M"][i] = mod
        if i + 3 < len(df["# ID"]) and str(df["Region"][i + 3]) == "lipid-modified cysteine":
            df["lipid_C"][i] = df["Region_L"][i + 3]

df = df[(df["N_pos1"]> 0)]
df=df.drop(columns=['Software', 'Region', 'Region_L', 'Region_R'])
df.to_excel(file_path3 + "/" + "region_output1_temp.xlsx", index=False)

#3 Reading the “prediction_results.txt” file
print("3 Reading the “prediction_results.txt” file")
df2 = pd.read_table(file_path1 + "/" +  "prediction_results.txt", skiprows=[0])
col_name = df2.columns.tolist()
add_list=["SeqID","SeqID_rev","sp_score"]
for i in range(len(add_list)):
    index = len(col_name) +i+1
    col_name.insert(index, add_list[i])
    df2=df2.reindex(columns=col_name)

for j in range(0, len(df2["Prediction"])):
    start3 = df2["# ID"][j].index("-")
    start5 = df2["# ID"][j].index("_P")
    df2["SeqID"][j] = df2["# ID"][j][:start3]
    df2["SeqID_rev"][j] = df2["# ID"][j][:start5]
    str_ = str(df2["CS Position"][j])
    if "Pr: " in str_:
        df2["sp_score"][j] = str_[str_.index("Pr: ") + 4:]
df2.to_excel(file_path3 + "/" +"prediction_results_temp.xlsx", index=False)

#4 document match
print("4 document match")
df3 = pd.read_excel(file_path3 + "/" + "region_output1_temp.xlsx")
df4 = pd.read_excel(file_path3 + "/" + "prediction_results_temp.xlsx")
for i in range(len(df3["SeqID_rev"])):
    if df3["lipid_C"][i] > 1:
        df3["pos_SP"][i] = df3["lipid_C"][i] + df3["pos_M"][i]
    else:
        df3["pos_SP"][i] = df3["C_pos2"][i]+ df3["pos_M"][i]
SeqID_ls=df4["SeqID_rev"]
score_ls=df4["sp_score"]
mat=dict(zip(SeqID_ls,score_ls))
for i in range(len(df3["SeqID_rev"])):
    df3["sp_score"][i]=mat[df3["SeqID_rev"][i]]
df3.to_excel(file_path3 + "/" + f"{specie}"+"_M_sp_temp.xlsx", index=False)

# Step 5 Extracting Signal Peptide Localization Information
print("5 Extracting Signal Peptide Localization Information")
file_in = open(file_path2 + "/"+ f"{specie}"+"_M.fasta", 'r')
out_file1 = open(file_path3 + "/"+ f"{specie}"+"_M_sp.fasta", 'w')
out_file2 = open(file_path3 + "/"+ f"{specie}"+"_M_mpep_temp.fasta", 'w')

fa_Con = file_in.read()
file_in.close()
every_fas = fa_Con.split(">")
title_con = []
seq_con = []
for i in every_fas:
    if i != "":
        start = i.index(" ")
        start1 = i.index("\n")
        title = i[:start]
        seq = i[start1:].replace("\n","")
        title_con.append(title)
        seq_con.append(seq)
fa = dict(zip(title_con, seq_con))
df5=pd.read_excel(file_path3 + "/"+ f"{specie}"+"_M_sp_temp.xlsx")
for j in range (len(df5["SeqID_rev"])):
    title = df5["SeqID_rev"][j]
    seq_con =fa[f"{title}"]
    N_pos1 = df5["N_pos1"][j]
    N_pos2 = df5["N_pos2"][j]
    H_pos1 = df5["H_pos1"][j]
    H_pos2 = df5["H_pos2"][j]
    C_pos1 = df5["C_pos1"][j]
    C_pos2 = df5["C_pos2"][j]
    lipd_C = df5["lipid_C"][j]
    df5["Seq_N"][j] = seq_con[int(N_pos1 - 1):int(N_pos2)]
    df5["Seq_H"][j] = seq_con[int(H_pos1 - 1):int(H_pos2)]
    df5["Seq_C"][j] = seq_con[int(C_pos1 - 1):int(C_pos2)]
    if lipd_C>0:
        df5["Seq_lipd"][j] = seq_con[int(C_pos2 - 1):int(lipd_C)]
    pos_sp=df5["pos_SP"][j]-df5["pos_M"][j]
    Sq_=fa[f"{title}"]
    Sq_sp=Sq_[pos_sp:]
    out_file1.write(">" + title+" "+ f"{specie}"+" " +str(len(Sq_))+"\n" +Sq_+"\n")
    out_file2.write(">" + title + " " + f"{specie}"+ " " + str(len(Sq_sp)) + "\n"+Sq_sp+"\n")
df5.to_excel(file_path3 + "/" + f"{specie}"+"_M_sp2_temp.xlsx", index=False)
file_in=open(file_path3 + "/"+ f"{specie}"+"_M_mpep_temp.fasta", 'r')
file_out=open(file_path3 + "/"+ f"{specie}"+"_M_mpep.fasta", 'w')
pro_der(file_in,file_out)

#6 Compute cysteine framework and output fasta document
print("6 Compute cysteine framework and output fasta document")
spscore_limit=0.9
df6 = pd.read_excel(file_path3 + "/"+ f"{specie}"+"_M_sp2_temp.xlsx")
file_in1 = open(file_path3 + "/"+ f"{specie}"+"_M_mpep.fasta", 'r')
fa_con1 = file_in1.read()
file_in1.close()
every_fas1 = fa_con1.split(">")
title_con1 = []
seq_con1 = []
for i in every_fas1:
    if i != "":
        start = i.index(" ")
        title = i[:start]
        start2 = i.index("\n")
        seq = i[start2:].replace('\n', '')
        title_con1.append(title)
        seq_con1.append(seq)
fa1 = dict(zip(title_con1, seq_con1))
fa1_index=list(fa1.keys())
fa1_index_old = [i[0:-2]for i in fa1_index]

df6["Sq_Len"] = ""
df6["C_num"] = ""
df6["C_Pos"] = ""
df6["C_frame"] = ""
for j in range(len(df6["SeqID_rev"])):
    t1= df6["SeqID_rev"][j]
    if t1 in fa1:
        sq1_ =fa1[t1]
        df6["Sq_Len"][j] = len(sq1_)
        df6["C_num"][j]= sq1_.count(seq_screen)
        st_1 = sq1_.find(seq_screen)
        a = []
        b = []
        c = 0
        if st_1 != -1 and st_1 != len(sq1_) - len(seq_screen):
            a = [m.start() for m in re.finditer(seq_screen, sq1_)]
            k = 0
            while k < len(a) - 1:
                c = a[k + 1] - a[k]
                b.append(c)
                k += 1
        df6["C_Pos"][j] = str(a)
        df6["C_frame"][j] = str(b)
    if t1 not in fa1:
        if t1[0:-2] in fa1_index_old:
            indices= fa1_index_old.index(t1[0:-2])
            # indices = [i for i, x in enumerate(fa1_index_old) if x == t1[0:-2]]
            t1_rev=fa1_index[indices]
            sq1_ = fa1[t1_rev]
            df6["Sq_Len"][j] = len(sq1_)
            df6["C_num"][j] = sq1_.count(seq_screen)
            st_1 = sq1_.find(seq_screen)
            a = []
            b = []
            c = 0
            if st_1 != -1 and st_1 != len(sq1_) - len(seq_screen):
                a = [m.start() for m in re.finditer(seq_screen, sq1_)]
                k = 0
                while k < len(a) - 1:
                    c = a[k + 1] - a[k]
                    b.append(c)
                    k += 1
            df6["C_Pos"][j] = str(a)
            df6["C_frame"][j] = str(b)

df6 = df6[(spscore_limit < df6["sp_score"])&(len(df6["Sq_Len"])>0)]
df6.to_excel(file_path3 + "/"+ f"{specie}"+"_M_sp.xlsx", index=False)

df7 = pd.read_excel(file_path3 + "/"+ f"{specie}"+"_M_sp.xlsx")
df7 = df7[(1 < df7["C_num"])]
df7.to_excel(file_path3 + "/"+ f"{specie}"+"_M_sp_2C.xlsx", index=False)
df8 = pd.read_excel(file_path3 + "/"+ f"{specie}"+"_M_sp_2C.xlsx")
out_file = open(file_path3 + "/"+ f"{specie}"+"_M_sp_2C_mpep.fasta", 'w')
for j in range(len(df8["SeqID_rev"])):
    t1= df8["SeqID_rev"][j]
    if t1 in fa1:
        sq1_ =fa1[t1]
        out_file.write(">" + t1 + " " + f"{specie}" + " " + str(len(sq1_)) + "\n" + sq1_ + "\n" )

# step7: Delete excess files
print("Step7: Delete excess files")
fname = os.listdir(file_path3)
fname = [i for i in fname if '_temp' in i]
for i in fname:
    os.remove(file_path3 + "/" + i)

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
