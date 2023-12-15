"""
Purpose:
1 Correcting mgf files within the range of -2, -1, and 0.
2 Extracting MS1 information from mgf files before and after correction.
3 Merging all Excel files.
"""
import csv
import re
import time
import os
import pandas as pd

start_time = time.time()
def output_list(list_tmp, m_origin,Title_origin):
    """
    object: output the mgf documents
    :param list_tmp: mgf document list
    :param m_origin: modified m/z
    :param Title_origin: title
    """
    for x in list_tmp:
        first_char = x[0][0]
        if str.isalpha(first_char):
            Title_ = re.search(r'(?<=TITLE=).*', ''.join(x[0]))
            if Title_ is not None:
                print("TITLE=", Title_origin, " ", x[1], " ", x[2], " ", x[3]," ", x[4], sep='', file=f_output)
                continue

            m_ = re.search(r'(?<=PEPMASS=).*', ''.join(x[0]))
            if m_ is not None:
                m = float(m_.group(0))
                if m > 0:
                    if x[1]!=None:
                        print("PEPMASS=", m_origin, " ", x[1], sep='', file=f_output)
                    else:
                        print("PEPMASS=", m_origin, " ", "0", sep='', file=f_output)
                    continue
        for i in range(len(x)):
            print(x[i],end=' ',file=f_output)
        print("", end='\n', file=f_output)

def mgf_ms1(file_path,file_name):
    """
    object：Extracting the MS1 information
    :param file_path: file path
    :param file_name: file name list
    :return:
    """
    for file in file_name:
        start = file.index(".")
        file_rev = file[:start]
        f_output = open(file_path + '/' + file_rev + "_title.csv", 'w')
        print("TITLE", "scan", "RTINSECONDS", "PEPMASS", "CHARGE", "data_M", "RT", "sample", end='\n', file=f_output)
        MS2_count = 0
        with open(file_path + '/' + file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            # line_count = 0
            for row in csv_reader:
                if row[0] == "BEGIN":
                    list_tmp = []
                    MS2_count = MS2_count + 1
                    print("##########--------------##########")
                    print(row[0], "this is the", MS2_count, "MS2 fragment ions")
                    RT = 0
                    m = 0
                    Title = []
                    scan = 0
                    RTinseconds = 0
                    data_M = 0
                    sample = []

                list_tmp.append(row)
                first_char = row[0][0]
                if str.isalpha(first_char):
                    Title_ = re.search(r'(?<=TITLE=).*', ''.join(row[0]))
                    if Title_ is not None:
                        Title = Title_.group(0)
                        T_s = Title
                        print("Title is:", Title)
                        if "/" in T_s:
                            start1 = T_s.index("/")
                            data_M = T_s[start1 + 1:]
                        else:
                            data_M = 0
                        start2 = T_s.index(".")
                        sample = T_s[:start2]

                    scan_ = re.search(r'(?<=scan=)\w+', ''.join(row[-1]))
                    if scan_ is not None:
                        scan = scan_.group(0)
                        # print("scan is:", scan)

                    RT_ = re.search(r'(?<=RTINSECONDS=).*', ''.join(row[0]))
                    if RT_ is not None:
                        RTinseconds = float(RT_.group(0))
                        RT = float(RT_.group(0)) / 60
                        # print("保留时间RT是:", RT, "秒")

                    m_ = re.search(r'(?<=PEPMASS=).*', ''.join(row[0]))
                    if m_ is not None:
                        m = float(m_.group(0))
                        # print("质荷比（m）是:", '%.4f' % m)

                    c_ = re.search(r'(?<=CHARGE=).*', ''.join(row[0]))
                    if c_ is not None:
                        c = str(c_.group(0))
                        print(Title, scan, RTinseconds, str(m), str(c), data_M, RT, sample, end='\n', file=f_output)
        f_output.close()

#step1. Reading documents
print("step1. Reading documents")
file_path = r'C:\Users\bqxin\Desktop\SmallP\raw data' # Folder path containing the raw mass spectrometry data
file_path1 = r'C:\Users\bqxin\Desktop\SmallP\mgf_cal' # Folder path of the results file
result_title="DL-TQY-fraction_RAE"
file_name = os.listdir(file_path)
file_name = [i for i in file_name if '.mgf' in i]
file_name.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))

#step2. monoisotope mass calibration
print("step2. monoisotope mass calibration")
for file in file_name:
    start = file.index(".")
    file_rev = file[:start]
    f_output = open(file_path1+'/'+file_rev+"_rev.mgf", 'w')
    MS2_count = 0
    with open(file_path+'/'+file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            if row[0] == "BEGIN":
                list_tmp = []
                MS2_count = MS2_count + 1
                print("##########--------------##########")
                print(row[0], "this is the", MS2_count, "MS2 fragment ions")

                # Initialize all parameters to zero
                RT = 0  # retention time
                m = 0  # m/z
                c = 0  # Charge
                M = 0  # mass weight
                Title = []
                scan = 0
                frag_M = []
                m_origin = 0
            list_tmp.append(row)

            first_char = row[0][0]
            if str.isalpha(first_char):
                Title_ = re.search(r'(?<=TITLE=).*', ''.join(row[0]))
                if Title_ is not None:
                    Title = Title_.group(0)
                    # print("Title is:", Title)

                scan_ = re.search(r'(?<=scan=)\w+', ''.join(row[-1]))
                if scan_ is not None:
                    scan = scan_.group(0)
                    # print("scan is:", scan)

                RT_ = re.search(r'(?<=RTINSECONDS=).*', ''.join(row[0]))
                if RT_ is not None:
                    RT = float(RT_.group(0))
                    # print("保留时间RT是:", RT, "秒")

                m_ = re.search(r'(?<=PEPMASS=).*', ''.join(row[0]))
                if m_ is not None:
                    m = float(m_.group(0))
                    # print("质荷比（m）是:", '%.4f' %m )

                c_ = re.search(r'(?<=CHARGE=)\w+', ''.join(row[0]))
                if c_ is not None:
                    c = int(c_.group(0))
                    # print("质荷比（c）是:", '%.0f' %m )

            if str.isdigit(first_char):
                frag_M.append(row[0])

            if row[0] == "END":
                M=m*c-1.0073*c
                if M>=2000:
                    if M>=4000:
                        for j in range(-2, 1):
                            M_rev = m * c - 1.0073 * c + j * 1.0073
                            m_origin = (M_rev + c * 1.0073) / c
                            Title_origin = Title + "/" + str(j)
                            output_list(list_tmp, m_origin, Title_origin)
                    else:#4000>M>=2000:
                        for j in range(-1, 1):
                            M_rev = m * c - 1.0073 * c + j * 1.0073
                            m_origin = (M_rev + c * 1.0073) / c
                            Title_origin = Title + "/" + str(j)
                            output_list(list_tmp, m_origin, Title_origin)
                else:#M<=2000:
                    m_origin = m
                    Title_origin = Title + "/" + str(0)
                    output_list(list_tmp, m_origin, Title_origin)
    f_output.close()

#Step3 extraction of MS1 information from the MGF file
print("Step3 extraction of MS1 information from the MGF file")
file_name1 = os.listdir(file_path1)
file_name1 = [i for i in file_name1 if '_rev.mgf' in i]
file_name1.sort(key=lambda x: str(''.join(filter(str.isalpha, x))))
mgf_ms1(file_path1,file_name1)

#step4 combination of the csv documents
print("step4 combination of the csv documents")
data_list=[]
fname = os.listdir(file_path1)
fname = [i for i in fname if '.csv' in i]
file_name.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
for i in fname:
    data_list.append(
        pd.read_csv(file_path1 + "/" + i, sep=" "))
data_all = pd.concat(data_list)
data_all.to_excel(file_path1 + "/" + f"{result_title}.xlsx", index=False)

# step5: Delete excess files
print("Step5: Delete excess files")
fname = os.listdir(file_path1)
fname = [i for i in fname if '.csv' in i]
for i in fname:
    os.remove(file_path1 + "/" + i)

# Computational time consumption
end_time = time.time()
escape = (end_time - start_time)/60
print('The code run {:.0f}m {:.0f}s'.format(escape // 60, escape % 60))
print("程序结束")

#————————————————
"""
Qirui Bi
Shanghai
December 7, 2023
"""