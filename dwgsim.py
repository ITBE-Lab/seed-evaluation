import os
import subprocess
import random
import multiprocessing
from config import *

def create_illumina_reads_dwgsim(sequenced_genome_path, reads_folder, num_reads, name, read_length):
    print("\tdwgsim...")
    dwgsim_instances = []
    file_names1 = ""
    file_names2 = ""
    reads1 = reads_folder + name + ".bwa.read1.fastq.gz"
    reads2 = reads_folder + name + ".bwa.read2.fastq.gz"
    # we actually need to seed dwgsim otherwise each instance will create the same data (it uses the current time)
    seed = random.randrange(0, 2**6)
    num_instances = 1#32
    for idx in range(num_instances):
        command = [dwgsim_str, "-1", str(read_length), "-2", str(read_length), "-N",
                    str(int(num_reads/num_instances)),  "-o", "1", "-z", str(seed), sequenced_genome_path,
                    reads_folder + "part_" + str(idx) + name]
        dwgsim_instances.append(subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        seed += 42 # increment the seed so that we have some different number every time
        file_names1 += reads_folder + "part_" + str(idx) + name + ".bwa.read1.fastq.gz "
        file_names2 += reads_folder + "part_" + str(idx) + name + ".bwa.read2.fastq.gz "
    for proc in dwgsim_instances:
        proc.wait()

    print("\tcat...")
    os.system("zcat " + file_names1 + " > " + reads1)
    os.system("zcat " + file_names2 + " > " + reads2)

def create_reads_survivor(sequenced_genome_path, reads_folder, num_reads, name, error_profile):
    print("\tsurvivor...")
    reads1 = reads_folder + name + ".fasta"

    command = survivor_str + sequenced_genome_path + " " + error_profile + " " + str(num_reads) + " " \
              + reads1 + " >/dev/null 2>&1"
    p = multiprocessing.Process(target=os.system, name="survivor command", args=(command,))
    p.daemon = True # kill p when this process dies...
    p.start()

    p.join(30)
    if p.is_alive():
        print("WARNING: survivor is running for 30 secs. It seems to be stuck; killing it now.")
        p.terminate()
        p.join()
