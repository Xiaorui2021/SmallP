"""
object：
1. Alignment RA and RAE group peptides with Nat group peptides.
2 Matching protein sequence information with Nat group peptides.
3 Matching signal peptide information with Nat group peptides.
"""
import pandas as pd
import os
import re
import time
import winsound

def Seq_match(df,df1,pos_1,pos_2,pos_3,pos_4,pos_5,pos_6,pos_7,pos_8,pos_9, pos_10, pos_11, pos_12, pos_13, pos_14, pos_15):
    num1 =0
    for i in range(len(df["_NO."])):
        seq_=df["_Sequence"][i]
        samp1=df["_Sample"][i]
        PTM_U1=df["_PTM-Unimod"][i]
        ls_seq1 = []
        ls_samp=[]
        ls_=[]
        for j in range (len(df1["_NO."])):
            seq_1 = df1["_Sequence"][j]
            samp2=df1["_Sample"][j]
            PTM_U2 = df1["_PTM-Unimod"][j]
            if samp1==samp2 and seq_.find(seq_1)!= -1 and len(PTM_U1)==0:
                ls_seq1.append (seq_1)
                ls_samp.append(samp2)
                ls_.append(j)
            PTM_U = set(PTM_U1) & set(PTM_U2)
            if samp1==samp2 and seq_.find(seq_1) != -1 and len(PTM_U1) > 0 and len(PTM_U)>0:
                ls_seq1.append(seq_1)
                ls_samp.append(samp2)
                ls_.append(j)
        if len(ls_seq1)>0:
            seq_max1 = max(ls_seq1, key=len, default='')
            k=[i for i, x in enumerate(ls_seq1) if x == max(ls_seq1)]
            for m in k:
                if ls_samp[m]==samp1:
                    a=ls_[m]
                    df[pos_1][i]=seq_max1
                    df[pos_2][i]=df1["_m/z"][a]
                    df[pos_3][i] = df1["_RT"][a]
                    df[pos_4][i] = df1["_Mass"][a]
                    df[pos_5][i] = df1["_Score"][a]
                    df[pos_6][i] = df1["_Ppm"][a]
                    df[pos_7][i] = df1["_length"][a]
                    df[pos_8][i] = df1["_Sample"][a]
                    df[pos_9][i] = df1["_Accession"][a]
                    df[pos_10][i] = df1["_PTM-Unimod"][a]
                    df[pos_11][i] = df1["_Scan"][a]
                    df[pos_12][i] = df1["_DetaM"][a]
                    df[pos_13][i] = df1["_Software"][a]
                    df[pos_14][i] = set(df["_PTM-Unimod"][i]) & set(df1["_PTM-Unimod"][a])
                    df[pos_15][i] = 100 * len(df1["_Sequence"][a]) / len(df["_Sequence"][i])

def pro_match(df,f_sp,n,dic1):
    for j in range(0,len(df["_NO."])):
        Seq_=df["_Sequence"][j]
        f_sp.seek(0)
        fa_Con = f_sp.read()
        every_fas = fa_Con.split(">")
        l_title=[]
        l_seq=[]
        l_socre = []
        for i in every_fas:
            if i != "":
                start = i.index("\n")
                start2 = i.index(" ")
                seq_con = i[start:].replace("\n","")
                seq_title = i[:start2]
                str_ = seq_con.find(Seq_)
                if str_ != -1:
                    l_title.append(seq_title)
                    l_seq.append(seq_con)
                    l_socre.append(str(round(dic1[seq_title],4)))
        p=0
        sp_max = ""
        if len(l_socre) >0:
            p = l_socre.index(max(l_socre))
            df["SeqID_rev"][j]=l_title[p]
            sp_max=l_seq[p]
        df["Sq_Len"][j] = len(sp_max)
        str_1 = sp_max.find(Seq_)
        df["Sq_pos_l"][j] = str_1+1
        df["Sq_pos_r"][j] = str_1+len(Seq_)
        N_left_pos = str_1 - n
        N_right_pos = str_1 + n
        C_pos = str_1 + len(Seq_)
        C_left_pos = C_pos - n
        C_right_pos = C_pos + n
        C_end = len(sp_max)
        if str_1 > n - 1 and str_1 + len(Seq_) + n <= len(seq_con):
            df["N_left"][j] = sp_max[N_left_pos:str_1]
            df["N_right"][j] = sp_max[str_1:N_right_pos]
            df["C_left"][j] = sp_max[C_left_pos:C_pos]
            df["C_right"][j] = sp_max[C_pos:C_right_pos]
        if str_1 <= n - 1:
            df["N_left"][j] = sp_max[start:str_1]
            df["N_right"][j] = sp_max[str_1:N_right_pos]
            df["C_left"][j] = sp_max[C_left_pos:C_pos]
            df["C_right"][j] = sp_max[C_pos:C_right_pos]
        if str_1 > n - 1 and str_1 + len(Seq_) + n > len(seq_con):
            df["N_left"][j] = sp_max[N_left_pos:str_1]
            df["N_right"][j] = sp_max[str_1:N_right_pos]
            df["C_left"][j] = sp_max[C_left_pos:C_pos]
            df["C_right"][j] = sp_max[C_pos:C_end]

#Step 1: Import data
print("第一步，导入数据")
start_time= time.time()

file_path= r'C:\Users\bqxin\Desktop\SmallP\Peaks_cal\cal'
file_path1= r'C:\Users\bqxin\Desktop\SmallP\raw data'
f_sp = open(file_path1 + "/" + "Pheretima aspergillum_M_sp_2C_mpep.fasta", 'r')
df_sp=pd.read_excel(file_path1 + "/" + "Pheretima aspergillum_M_sp.xlsx", sheet_name=0)
seqIDrev_l=df_sp["SeqID_rev"]
spscore_l=df_sp["sp_score"]
dic1=dict(zip(seqIDrev_l,spscore_l))
tl="Pheretima aspergillum-sp-Nat_SPIDER_identification"

n=5
ppm_low=-5
ppm_high=10
df = pd.read_excel(file_path + "/" + tl+".xlsx", sheet_name=0)

# Step2: Alignment RA and RAE group peptides with Nat group peptides.
print("Step2: Alignment RA and RAE group peptides with Nat group peptides.")
file_name = os.listdir(file_path)
file_name = [i for i in file_name if '_identification.xlsx' in i]
file_name.sort(key=lambda x: str(''.join(filter(str.isalpha, x))))
for file in file_name:
    print(file)
    species=(re.split(r'[-_]', file))[0]
    db = (re.split(r'[-_]', file))[1]
    type= (re.split(r'[-_]', file))[2]
    software = (re.split(r'[-_]', file))[3]
    document = species + "-" + db + "-" + type + "_" + software
    if type !="Nat":
        print(type)
        col_name = df.columns.tolist()
        pos_1 = f"{type}_seq"
        pos_2 = f"{type}_m/z"
        pos_3 = f"{type}_RT"
        pos_4 = f"{type}-Mass"
        pos_5 = f"{type}-Score"
        pos_6 = f"{type}-ppm"
        pos_7 = f"{type}-Len"
        pos_8 = f"{type}-sample"
        pos_9 = f"{type}-Accession"
        pos_10 = f"{type}-PTM"
        pos_11 = f"{type}-Scan"
        pos_12 = f"{type}-deta_M"
        pos_13 = f"{type}-software"
        pos_14 = f"{type}-PTM_U"
        pos_15 = f"{type}-rate"
        add_list = [pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, pos_7, pos_8, pos_9, pos_10, pos_11, pos_12, pos_13,
                    pos_14, pos_15]
        for i in range(len(add_list)):
            index = len(col_name) + i + 1
            col_name.insert(index, add_list[i])
            df = df.reindex(columns=col_name)
        df1 = pd.read_excel(file_path + "/" + document+"_identification.xlsx", sheet_name=0)
        Seq_match(df, df1, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, pos_7, pos_8, pos_9, pos_10, pos_11, pos_12,
                  pos_13, pos_14, pos_15)

df.to_excel(file_path + "/" + "RA-RAE matched_temp.xlsx", index=False)
df2 = df[(df["_Ppm"] >ppm_low) & (df["_Ppm"] <ppm_high)]
df2 = df2[(df2["RA-Len"] >0) | (df2["RAE-Len"] >0)]
df2.to_excel(file_path + "/" + "refined peptides_temp.xlsx", index=False)

# step3 Matching protein sequence information with Nat group peptides
print(" step3 Matching protein sequence information with Nat group peptides")
df = pd.read_excel(file_path + "/" + "refined peptides_temp.xlsx", sheet_name=0)
col_name = df.columns.tolist()
pos_1 = "SeqID"
pos_2 = "SeqID_rev"
pos_3 = "Sq_Len"
pos_4 = "N_left"
pos_5 = "N_right"
pos_6 = "C_left"
pos_7 = "C_right"
pos_8 = "Seq_N"
pos_9 = "Seq_H"
pos_10 = "Seq_C"
pos_11 = "Seq_lipd"
pos_12= "SP_pos_l"
pos_13= "SP_pos_r"
pos_14="Sq_pos_l"
pos_15="Sq_pos_r"
pos_16="M_pos"
pos_17="deta_pos"
pos_18="sp_score"
pos_21="Pro_Cnum"
pos_22="Pro_Cframe"
add_list = [pos_1, pos_2,pos_3, pos_4,pos_5, pos_6,pos_7, pos_8,pos_9, pos_10, pos_11,
            pos_12, pos_13, pos_14, pos_15, pos_16, pos_17, pos_18, pos_21, pos_22]
for i in range(len(add_list)):
    index = len(col_name) + i + 1
    col_name.insert(index, add_list[i])
    df = df.reindex(columns=col_name)
pro_match(df,f_sp,n,dic1)
df.to_excel(file_path + "/" + "RA-RAE-NAE-sp matched_temp.xlsx", index=False)

#step4 Matching signal peptide information with Nat group peptides.
print("step4 Matching signal peptide information with Nat group peptides")
df1=pd.read_excel(file_path + "/" + "RA-RAE-NAE-sp matched_temp.xlsx", sheet_name=0)
for i in range (len(df1["_Sequence"])):
    for j in range (len(df_sp["SeqID_rev"])):
        if df1["SeqID_rev"][i]==df_sp["SeqID_rev"][j]:
            df1['Seq_N'][i]=df_sp['Seq_N'][j]
            df1['Seq_H'][i] = df_sp['Seq_H'][j]
            df1['Seq_C'][i] = df_sp['Seq_C'][j]
            df1['Seq_lipd'][i] =df_sp['Seq_lipd'][j]
            df1['SP_pos_l'][i] = 1
            df1['SP_pos_r'][i] = df_sp['pos_SP'][j]-df_sp['pos_M'][j]
            df1['M_pos'][i] = df_sp['pos_M'][j]
            df1['Pro_Cnum'][i] = df_sp['C_num'][j]
            df1['Pro_Cframe'][i] = df_sp['C_frame'][j]
            df1['deta_pos'][i]=df1['Sq_pos_l'][i]-df1['SP_pos_r'][i]-1
            df1['sp_score'][i] = df_sp['sp_score'][j]
df1.to_excel(file_path + "/" + tl+"_peptides.xlsx", index=False)

# step5: Delete excess files
print("Step5: Delete excess files")
fname = os.listdir(file_path)
fname = [i for i in fname if '_temp' in i]
for i in fname:
    os.remove(file_path + "/" + i)

# Computational time consumption
end_time = time.time()
escape = end_time - start_time
print('The code run {:.0f}m {:.0f}s'.format(escape // 60, escape % 60))
print("程序结束")
#————————————————
"""
Qirui Bi
Shanghai
December 7, 2023
"""

winsound.PlaySound('SystemHand', winsound.SND_ALIAS)