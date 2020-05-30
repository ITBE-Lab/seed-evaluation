import os
import subprocess
import random
import multiprocessing
from config import *
from bokeh.plotting import figure, show

def create_illumina_reads_dwgsim(sequenced_genome_path, reads_folder, num_reads, name, read_length, error_factor=1):
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
                    str(int(num_reads/num_instances)), "-r", str(0.0010*error_factor), "-o", "1", "-z",
                    str(seed), sequenced_genome_path,
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

def gen_survivor_error_profile(out_file_name, in_file_name=survivor_error_profile, p_mod=lambda x: x):
    with open(os.path.expanduser(in_file_name)) as in_file:
        lines = in_file.readlines()[1:]
        pos_list = [int(line.split("\t")[0]) for line in lines]
        p_stop = [float(line.split("\t")[1]) for line in lines]
        p_match = [float(line.split("\t")[2]) for line in lines]
        p_missmatch = [float(line.split("\t")[3]) for line in lines]
        p_ins = [float(line.split("\t")[4]) for line in lines]
        p_del = [float(line.split("\t")[5]) for line in lines]
    with open(out_file_name, "w") as out_file:
        out_file.write("Pos\tP(stop)\tP(match)\tp(mismatch)\tP(ins)]\tP(del)\n")
        for pos in pos_list:
            p_stop_val = p_stop[pos]
            p_missmatch_val = p_mod(p_missmatch[pos])
            p_match_val = 1 - p_mod(1 - p_match[pos])
            p_ins_val = p_mod(p_ins[pos])
            p_del_val = p_mod(p_del[pos])
            out_file.write(str(pos) + "\t" + str(p_stop_val) + "\t" + str(p_match_val) + "\t" + str(p_missmatch_val) +
                           "\t" + str(p_ins_val) + "\t" + str(p_del_val) + "\n")

def gen_survivor_error_profile_fac(out_file_name, fac=1, survivor_error_profile=survivor_error_profile):
    gen_survivor_error_profile(out_file_name + "_" + str(fac) + ".txt",
                               survivor_error_profile,
                               p_mod=lambda x: min(1,x*fac))

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

def plot_error_profile(in_file_name):
    plot = figure(title=in_file_name)
    with open(os.path.expanduser(in_file_name)) as in_file:
        x = []
        ys = {"stop":[], "match": [], "ins": [], "del": []}
        for line in in_file.readlines()[1:]:
            pos, stop, match, _, ins, del_ = [float(x) for x in line.split("\t")]
            x.append(pos)
            ys["stop"].append(stop)
            ys["match"].append(match)
            ys["ins"].append(ins)
            ys["del"].append(del_)
        plot.line(x=x, y=ys["stop"], line_color="black", legend_label="stop")
        plot.line(x=x, y=ys["match"], line_color="red", legend_label="match")
        plot.line(x=x, y=ys["ins"], line_color="yellow", legend_label="ins")
        plot.line(x=x, y=ys["del"], line_color="green", legend_label="del")

    show(plot)

if __name__ == "__main__":
    #plot_error_profile(survivor_error_profile)
    gen_survivor_error_profile_fac("test_error_profile", 1 )
    plot_error_profile("test_error_profile_1.txt")