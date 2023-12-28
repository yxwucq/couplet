from typing import List

import re
import subprocess
import os
import itertools
import gzip
import sys
import time

# modified from https://github.com/yupenghe/methylpy/blob/methylpy/methylpy/utilities.py
def split_fastq_file(num_chunks: int,
                     input_file: str) -> int:
    """
    This function mimics the unix split utility. 
    """
    
    file_handles = {}
    
    for index in range(0, num_chunks):
        if input_file.endswith('_R1.fq.gz'):
            outname = input_file.replace('_R1.fq.gz', f"_SPLIT_{index}_R1.fq.gz")
        elif input_file.endswith('_R2.fq.gz'):
            outname = input_file.replace('_R2.fq.gz', f"_SPLIT_{index}_R2.fq.gz") 
        file_handles[index] = gzip.open(outname,'wt')
        
    cycle = itertools.cycle(list(range(0,num_chunks)))
    total_reads = 0
    
    if input_file[-3:] == ".gz":
        f = gzip.open(input_file,'rt')
    else:
        f = open(input_file,'r')

    while True:
        current_file = next(cycle)
        lines = [f.readline() for _ in range(4)]
        if not lines[0]:
            break
        total_reads += 1
        file_handles[current_file].writelines(('').join(lines))
    f.close()

    for index in range(0,num_chunks):
        file_handles[index].close()

    return(total_reads)

def shell_split_fastq_file(input_file: str, 
                           chunks_size: int=1000000):
    
    # start_time = time.time()
    if input_file.endswith('_R1.fq.gz'):
        command = f"zcat {input_file} | parallel --pipe -N{chunks_size} 'gzip > {input_file.replace('_R1.fq.gz', '_SPLIT_')}{{#}}_R1.fq.gz'"
    elif input_file.endswith('_R2.fq.gz'):
        command = f"zcat {input_file} | parallel --pipe -N{chunks_size} 'gzip > {input_file.replace('_R2.fq.gz', '_SPLIT_')}{{#}}_R2.fq.gz'"
    else:
        command = f"zcat {input_file} | parallel --pipe -N{chunks_size} 'gzip > {input_file.replace('.fq.gz', '_SPLIT_')}{{#}}.fq.gz'"
    # print(command)
    subprocess.run(command, shell=True, check=True)
    # end_time = time.time()
    # print("Splitting finished in {} seconds.".format(end_time - start_time))
    if input_file.endswith('_R1.fq.gz'):
        split_file_list = [x for x in os.listdir(os.path.dirname(input_file)) if re.match(r'^.*_SPLIT_[0-9]+_R1.fq.gz$', x) and x.startswith(os.path.basename(input_file).replace('_R1.fq.gz', '_SPLIT_'))]
    elif input_file.endswith('_R2.fq.gz'):
        split_file_list = [x for x in os.listdir(os.path.dirname(input_file)) if re.match(r'^.*_SPLIT_[0-9]+_R2.fq.gz$', x) and x.startswith(os.path.basename(input_file).replace('_R2.fq.gz', '_SPLIT_'))]
    else:
        split_file_list = [x for x in os.listdir(os.path.dirname(input_file)) if re.match(r'^.*_SPLIT_[0-9]+.fq.gz$', x) and x.startswith(os.path.basename(input_file).replace('.fq.gz', '_SPLIT_'))]
    split_file_list = [os.path.join(os.path.dirname(input_file), x) for x in split_file_list]
    
    return split_file_list

def merge_fastq_files(input_files: List[str],
                      output_file: str) -> None:
    """
    This function mimics the unix cat utility.
    """
    if output_file[-3:] == ".gz":
        f = gzip.open(output_file,'wt')
    else:
        f = open(output_file,'w')
    
    for input_file in input_files:
        if input_file[-3:] == ".gz":
            g = gzip.open(input_file,'rt')
        else:
            g = open(input_file,'r')
        f.writelines(g.readlines())
        g.close()
        os.remove(input_file)
    f.close()
    
def shell_merge_fastq_files(input_files: List[str],
                            output_file: str):
    # start_time = time.time()
    command = f"cat {' '.join(input_files)} > {output_file}"
    # print(command)
    subprocess.run(command, shell=True, check=True)
    # end_time = time.time()
    # print("Merging finished in {} seconds.".format(end_time - start_time))