from MA import *
import math
import time
import dwgsim
from bokeh.plotting import figure, show, reset_output
from bokeh.palettes import Category20, Category10
from bokeh.io import export_png, export_svgs
import os
from bokeh_style_helper import *
from config import *
from runBackwardMEMOhlebusch import *

str_mm = str(get_mmi_parameter_set().by_name("Minimizer Window Size").get()) + "," + \
         str(get_mmi_parameter_set().by_name("Minimizer Size").get())

local_max_ambiguity_fmd = max_ambiguity_fmd

def create_genome_of_size(ref_pack, size, backwards=False):
    print("size:", size)
    contig_id = 0
    if backwards:
        contig_id = len(ref_pack.contigNames()) - 1
    pack_list = []
    while size > ref_pack.length_of_sequence_with_id(contig_id):
        #print(contig_id, ref_pack.start_of_sequence_id(contig_id), ref_pack.length_of_sequence_with_id(contig_id))
        nuc_seq = ref_pack.extract_from_to(ref_pack.start_of_sequence_id(contig_id),
                            ref_pack.start_of_sequence_id(contig_id) + ref_pack.length_of_sequence_with_id(contig_id))
        pack_list.append((ref_pack.name_of_sequence_id(contig_id), nuc_seq))
        size -= ref_pack.length_of_sequence_with_id(contig_id)
        if backwards:
            contig_id -= 1
        else:
            contig_id += 1
        del nuc_seq

    # center the extracted sequence in the contig
    start = (ref_pack.length_of_sequence_with_id(contig_id) - size) // 2 + ref_pack.start_of_sequence_id(contig_id)
    nuc_seq = ref_pack.extract_from_to(start, start + size)
    pack_list.append((ref_pack.name_of_sequence_id(contig_id), nuc_seq))

    if backwards:
        pack_list.reverse()

    new_pack = Pack()
    for a, b in pack_list:
        new_pack.append(a, "no description", b)

    return new_pack


def log_range(start=start_size, stop=stop_size, num_steps=num_steps):
    start = math.log(start)
    stop = math.log(stop)
    for x in range(0, num_steps, 1):
        yield round(math.exp(start + (stop-start)*x/(num_steps-1)))

def linear_range(start=start_size, stop=stop_size, num_steps=num_steps):
    if num_steps == 1:
        yield float(start)
    else:
        for x in range(0, num_steps, 1):
            yield start + (stop - start) * x / (num_steps - 1)

def set_up_folders():
    if not os.path.exists(prefix + "svg/"):
        os.mkdir(prefix + "svg/")
    if not os.path.exists(prefix + "genomes/"):
        os.mkdir(prefix + "genomes/")
    if not os.path.exists(prefix + "reads/"):
        os.mkdir(prefix + "reads/")
    if not os.path.exists(prefix + "profiles/"):
        os.mkdir(prefix + "profiles/")

def generate_genomes(time_steps, ref_pack, out_prefix, backwards=False):
    set_up_folders()
    print("generating genomes...")
    for idx, x in enumerate(time_steps()):
        print(idx, "...")
        ref_slice = create_genome_of_size(ref_pack, x, backwards)
        ref_slice.store(out_prefix + "/slice_" + str(x))
        with open(out_prefix + "/slice_" + str(x) + ".fasta", "w") as fasta_out:
            for name, sequence in zip(ref_slice.contigNames(), ref_slice.contigSeqs()):
                print("writing:", name, "len:", len(sequence))
                fasta_out.write(">")
                fasta_out.write(name)
                fasta_out.write("\n")
                idx = 0
                while idx < len(sequence):
                    fasta_out.write(sequence[idx:idx+50])
                    fasta_out.write("\n")
                    idx += 50
    print("done")

def generate_profiles(time_steps, survivor_error_profile=survivor_error_profile):
    set_up_folders()
    print("generating profiles...")
    for idx, x in enumerate(time_steps()):
        dwgsim.gen_survivor_error_profile_fac(prefix + "profiles/pacb", fac=x,
                                              survivor_error_profile=survivor_error_profile)
    print("done")

def generate_reads(time_steps, out_prefix, genome_prefix):
    print("generating reads...")
    for idx, x in enumerate(time_steps()):
        print(idx, "...")
        if x_axis_unit == "genome_section_size":
            dwgsim.create_illumina_reads_dwgsim(genome_prefix + "/slice_" + str(x) + ".fasta", out_prefix,
                                                num_illumina_reads, "slice_" + str(x), illumina_read_size)
            dwgsim.create_reads_survivor(genome_prefix + "/slice_" + str(x) + ".fasta", out_prefix, num_pacb_reads,
                                        "slice_" + str(x), survivor_error_profile)
        if x_axis_unit == "read_noise":
            if True:
                dwgsim.create_illumina_reads_dwgsim(reference_genome_fasta, out_prefix, num_illumina_reads,
                                                    "noise_" + str(x), illumina_read_size, error_factor=x )
            if True:
                dwgsim.create_reads_survivor(reference_genome_fasta, out_prefix, num_pacb_reads,
                                            "noise_" + str(x), prefix + "profiles/pacb_" + str(x) + ".txt")
    print("done")

def measure_time_section_size(caller, prefix, log_file_name, with_single=True, time_steps=log_range):
    print(log_file_name, ": measuring time...")
    with open(prefix + "illumina_" + log_file_name, "w") as log_file:
        log_file.write("index\tgenome size\truntime\n")
        for idx, x in enumerate(time_steps()):
            print(idx, "...")
            caller.prep(x, prefix)
            start = time.time()
            caller.run()
            end = time.time()
            caller.post()
            print("time required (paired):", end - start)
            log_file.write(str(idx) + "\t" + str(x) + "\t" + str(end - start) + "\n")
    if with_single:
        with open(prefix + "pacb_" + log_file_name, "w") as log_file:
            log_file.write("index\tgenome size\truntime\n")
            for idx, x in enumerate(time_steps()):
                print(idx, "...")
                caller.prep(x, prefix, paired=False)
                start = time.time()
                caller.run()
                end = time.time()
                caller.post()
                print("time required (pacbio):", end - start)
                log_file.write(str(idx) + "\t" + str(x) + "\t" + str(end - start) + "\n")

    print("done")

def measure_time_noise(caller, prefix, log_file_name, with_paired=True, with_single=True, time_steps=linear_range):
    ref_pack = Pack()
    ref_pack.load(reference_genome_path)
    p_m = get_mmi_parameter_set()
    mm_index = libMA.MinimizerIndex(p_m, prefix + "genomes/full.mmi")
    fm_index = FMIndex()
    fm_index.load(reference_genome_path)
    print(log_file_name, ": measuring time...")
    if with_paired:
        with open(prefix + "illumina_" + log_file_name, "w") as log_file_paired:
            log_file_paired.write("index\tgenome size\truntime\n")
            for idx, x in enumerate(time_steps()):
                print(idx, "...")

                reads_file_name = prefix + "reads/noise_" + str(x) + ".bwa.read1.fastq.gz"
                reads = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads_file_name)])) \
                            .cpp_module.read_all()

                caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                start = time.time()
                caller.run()
                end = time.time()
                caller.post()

                print("time required (paired):", end - start)
                log_file_paired.write(str(idx) + "\t" + str(x) + "\t" + str(end - start) + "\n")
    if with_single:
        with open(prefix + "pacb_" + log_file_name, "w") as log_file_single:
            log_file_single.write("index\tgenome size\truntime\n")
            for idx, x in enumerate(time_steps()):
                print(idx, "...")

                reads_file_name = prefix + "reads/noise_" + str(x) + ".fasta"
                reads = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads_file_name)])) \
                                .cpp_module.read_all()

                caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                start = time.time()
                caller.run()
                end = time.time()
                caller.post()

                print("time required (pacbio):", end - start)
                log_file_single.write(str(idx) + "\t" + str(x) + "\t" + str(end - start) + "\n")

    print("done")


def measure_time(caller, prefix, log_file_name, with_single=True):
    if x_axis_unit == "genome_section_size":
        measure_time_section_size(caller, prefix, log_file_name, with_single=with_single)
    if x_axis_unit == "read_noise":
        measure_time_noise(caller, prefix, log_file_name, with_paired=True, with_single=with_single)

def render_times(title, prefix, element_list, out_file, yAxisKey="runtime", y_axis_log=False, yAxisKey2=None,
                 yAxisKeyDivident=None, divide_y_by=1, y_range=None):
    xAxes = []
    yAxes = []
    legends = []
    colors = []
    for stack_list in element_list:
        currY = None
        for name, file_name, color in stack_list:
            legends.append(name)
            colors.append(color_scheme(color))
            with open(prefix + file_name, "r") as tsv_file:
                lines = [x[:-1].split("\t") for x in tsv_file.readlines()]
                if currY is None:
                    currY = [0]*(len(lines) - 1)
                header = {y:x for x,y in enumerate(lines[0])}
                xAxis_ = [ float(x[header["genome size"]]) for x in lines[1:] ]
                #assert xAxis is None or xAxis == xAxis_
                if yAxisKeyDivident is None:
                    yAxes_ = [ float(x[header[yAxisKey]])/divide_y_by + y for x,y in zip(lines[1:], currY) ]
                else:
                    yAxes_ = [ float(x[header[yAxisKey]])/(divide_y_by*float(x[header[yAxisKeyDivident]])) + y for x,y in zip(lines[1:], currY) ]
                xAxes.append(xAxis_)
                yAxes.append(yAxes_)
                currY = yAxes_
                if not yAxisKey2 is None:
                    yAxes_2 = [ float(x[header[yAxisKey2]])/divide_y_by for x in lines[1:] ]
                    yAxes.append(yAxes_2)
    reset_output()
    if x_axis_unit == "genome_section_size":
        if y_axis_log:
            plot = figure(title=title, x_axis_type="log", y_axis_type="log", y_range=y_range)
        else:
            plot = figure(title=title, x_axis_type="log", y_range=y_range)
    else:
        if y_axis_log:
            plot = figure(title=title, y_axis_type="log", y_range=y_range)
        else:
            plot = figure(title=title, y_range=y_range)
    step = 1 if yAxisKey2 is None else 2
    for yAxis, xAxis, legend, color in zip(yAxes[::step], xAxes, legends, colors):
        print(xAxis, yAxis)
        plot.line(x=xAxis, y=yAxis, legend_label=legend, line_color=color, line_width=point_to_px(3))
        plot.x(x=xAxis, y=yAxis, legend_label=legend, line_color=color, size=point_to_px(5), line_width=point_to_px(2))
        #plot.x(x=xAxis, y=yAxis, legend_label=legend, color=color)
    if not yAxisKey2 is None:
        for yAxis, xAxis, color in zip(yAxes[1::2], xAxes, colors):
            plot.line(x=xAxis, y=yAxis, line_color=color, line_dash=[1,1], line_width=point_to_px(3))
            #plot.cross(x=xAxis, y=yAxis, color=color)
    plot.legend.location = "top_left"
    if x_axis_unit == "genome_section_size":
        plot.xaxis.axis_label = "section length [nt]"
    if x_axis_unit == "read_noise":
        plot.xaxis.axis_label = "read noise"
    plot.yaxis.axis_label = yAxisKey
    style_plot(plot)
    show(plot)
    if save_plots:
        plot.output_backend = "svg"
        export_svgs(plot, filename=prefix + "svg/" + out_file + ".svg")


def get_read_positions(reads, pack):
    ret = []
    for read in reads:
        name = read.name.split("_")
        if name[0] == "rand" and name[1] == "0":
            ret.append((0, 0, True))
        else:
            start = pack.start_of_sequence(name[0]) + int(float(name[1])) - 1
            if len(name) > 3:
                strand = name[3] == "0"
            else:
                strand = name[2] == "+"
            ret.append((start, start+len(read), strand))
    return ret

## overlap between two intervals
def overlap(start_a, end_a, start_b, end_b):
    start = max(start_a, start_b)
    end = min(end_a, end_b)
    return max(0, end - start)


def plot_seeds(seeds, name, i_start, i_end):
    reset_output()
    plot = figure(title="seeds - " + name, plot_width=800, plot_height=800)
    xs = []
    ys = []
    size = i_end - i_start
    plot.rect(x=i_start + size/2, y=size/2, width=size, height=size, color="green")
    for seed in seeds:
        xs.append([seed.start_ref, seed.start_ref + seed.size])
        ys.append([seed.start, seed.start + seed.size])
    plot.multi_line(xs=xs, ys=ys)
    show(plot)

def get_seed_entropy(seeds_list, reads, pack, also_return_percent_covered=False, return_seeds_in_area=False,
                     name="unnamed"):
    if seeds_list is None:
        if also_return_percent_covered:
            return float("NaN"), 0
        if return_seeds_in_area:
            return 0, 0
        return 0
    read_pos = get_read_positions(reads, pack)
    # accumulated overlap with read region
    hits = 0
    # accumulated overlap with non read region
    num_nuc = 0
    # accumulated read size
    r_sum = 0
    # accumulated number of seeds hitting the read area
    seed_hits = 0

    num_seeds = 0
    for read in reads:
        r_sum += len(read)
    for (read_start, read_end, is_forward_strand), seeds in zip(read_pos, seeds_list):
        seeds.sort_by_ref_pos()
        max_r_pos = 0
        num_seeds += len(seeds)
        #plot_seeds(seeds, name, read_start, read_end)
        for seed in seeds:
            # this line makes sure we do not count overlaps with the same region twice (seeds are sorted by refpos)
            seed_start = max(seed.start_ref, max_r_pos)
            seed_end = seed.start_ref + seed.size

            max_r_pos = max(max_r_pos, seed_end)

            num_nuc += overlap(0, read_start, seed_start, seed_end)
            num_nuc += overlap(read_end, pack.unpacked_size(), seed_start, seed_end)
            if seed_end >= read_start and seed.start_ref <= read_end:
                seed_hits += 1
            hits += overlap(read_start, read_end, seed_start, seed_end)

    if also_return_percent_covered:
        if num_nuc == 0:
            return float("NaN"), 0
        return num_nuc / (pack.unpacked_size() * len(seeds_list)), hits / r_sum
    else:
        #if num_nuc + hits == 0:
        #    return 0
        #return hits / (num_nuc + hits)
        if num_seeds == 0:
            if return_seeds_in_area:
                return 0, 0
            return 0
        if return_seeds_in_area:
            return seed_hits / num_seeds, hits / num_seeds
        return hits / num_seeds

class CreateMinimizerIndex:
    def __init__(self):
        self.genome_prefix = ""
        self.contigs = None
        self.contig_names = None
        self.index = None
        self.p_m = None

    def prep(self, x, prefix, paired=True):
        self.genome_prefix = prefix + "genomes/slice_" + str(x)
        pack = Pack()
        pack.load(self.genome_prefix)
        self.contigs = pack.contigSeqs()
        self.contig_names = pack.contigNames()
        self.p_m = get_mmi_parameter_set()

    def run(self):
        self.index = libMA.MinimizerIndex(self.p_m, self.contigs, self.contig_names)
        self.index.set_mid_occ(local_max_ambiguity_fmd)

    def post(self):
        self.index.dump(self.genome_prefix + ".mmi")

class ComputeMinimizers:
    def __init__(self):
        self.index = None
        self.reads = None
        self.str_reads = None
        self.seeds = None
        self.ref_pack = None

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = get_mmi_parameter_set()
        self.index = libMA.MinimizerIndex(p_m, prefix + "genomes/slice_" + str(x) + ".mmi")
        self.index.set_mid_occ(local_max_ambiguity_fmd)
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        self.reads = reader.cpp_module.read_all()
        self.str_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            self.str_reads.append(str(nuc_seq))

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        self.ref_pack = ref_pack
        self.index = mm_index
        self.reads = reads
        self.str_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            self.str_reads.append(str(nuc_seq))

    def run(self):
        self.seeds = self.index.seed(self.str_reads, self.ref_pack)
        return self.seeds

    def post(self):
        pass
        #for seed, read in zip(self.seeds, self.reads):
        #    seed.confirm_seed_positions(read, self.ref_pack, False)

class LumpMinimizers:
    def __init__(self):
        self.seeds = None
        self.lumper = None
        self.nuc_seq_reads = None
        self.ref_pack = None
        self.lumped_seeds = None

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = get_mmi_parameter_set()
        index = libMA.MinimizerIndex(p_m, prefix + "genomes/slice_" + str(x) + ".mmi")
        self.index.set_mid_occ(local_max_ambiguity_fmd)
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        self.nuc_seq_reads = reader.cpp_module.read_all()
        reads = libMA.StringVector()
        for nuc_seq in self.nuc_seq_reads:
            reads.append(str(nuc_seq))
        self.seeds = index.seed(reads, self.ref_pack)
        self.lumper = libMA.SeedLumping(p_m)

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        p_m = get_mmi_parameter_set()
        self.ref_pack = ref_pack
        self.nuc_seq_reads = reads
        reads = libMA.StringVector()
        for nuc_seq in self.nuc_seq_reads:
            reads.append(str(nuc_seq))
        self.seeds = mm_index.seed(reads, self.ref_pack)
        self.lumper = libMA.SeedLumping(p_m)

    def run(self):
        self.lumped_seeds = self.lumper.cpp_module.lump(self.seeds, self.nuc_seq_reads, self.ref_pack)

    def post(self):
        x = 0
        for seeds in self.lumped_seeds:
            x += len(seeds)
        print("num seeds:", x)
        pass

class ExtendMinimizers:
    def __init__(self, min_seed_length=None):
        self.seeds = None
        self.ref_pack = None
        self.reads = None
        self.extender = None
        self.smems = None
        self.min_seed_length = min_seed_length

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = get_mmi_parameter_set()
        index = libMA.MinimizerIndex(p_m, prefix + "genomes/slice_" + str(x) + ".mmi")
        self.index.set_mid_occ(local_max_ambiguity_fmd)
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        self.reads = reader.cpp_module.read_all()
        std_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            std_reads.append(str(nuc_seq))
        minimizer_seeds = index.seed(std_reads, self.ref_pack)
        lumper = libMA.SeedLumping(p_m)
        self.seeds = lumper.cpp_module.lump(minimizer_seeds, self.reads, self.ref_pack)
        del minimizer_seeds
        self.extender = libMA.SeedExtender(p_m)

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        p_m = get_mmi_parameter_set()
        self.ref_pack = ref_pack
        self.reads = reads
        std_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            std_reads.append(str(nuc_seq))
        minimizer_seeds = mm_index.seed(std_reads, self.ref_pack)
        lumper = libMA.SeedLumping(p_m)
        self.seeds = lumper.cpp_module.lump(minimizer_seeds, self.reads, self.ref_pack)
        del minimizer_seeds
        self.extender = libMA.SeedExtender(p_m)

    def run(self):
        self.smems = self.extender.cpp_module.extend(self.seeds, self.reads, self.ref_pack)
        if not self.min_seed_length is None:
            return libMA.MinLength(get_mmi_parameter_set(), self.min_seed_length).cpp_module.filter(self.smems)
        return self.smems

    def post(self):
        pass
        #for seed, read in zip(self.smems, self.queries):
        #    seed.confirm_seed_positions(read, self.ref_pack, True)

class ExtendThenSortMinimizers:
    def __init__(self):
        self.seeds = None
        self.ref_pack = None
        self.reads = None
        self.extender = None
        self.smems = None

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = get_mmi_parameter_set()
        index = libMA.MinimizerIndex(p_m, prefix + "genomes/slice_" + str(x) + ".mmi")
        self.index.set_mid_occ(local_max_ambiguity_fmd)
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        self.reads = reader.cpp_module.read_all()
        std_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            std_reads.append(str(nuc_seq))
        self.minimizer_seeds = index.seed(std_reads, self.ref_pack)
        self.extender = libMA.SeedExtender(p_m)
        self.filter = libMA.SortRemoveDuplicates(p_m)

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        p_m = get_mmi_parameter_set()
        self.ref_pack = ref_pack
        self.reads = reads
        std_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            std_reads.append(str(nuc_seq))
        self.minimizer_seeds = mm_index.seed(std_reads, self.ref_pack)
        self.extender = libMA.SeedExtender(p_m)
        self.filter = libMA.SortRemoveDuplicates(p_m)

    def run(self):
        smems_w_dups = self.extender.cpp_module.extend(self.minimizer_seeds, self.reads, self.ref_pack)
        self.smems = self.filter.cpp_module.filter(smems_w_dups)
        return self.smems

    def post(self):
        pass

class MinimizerToSmem:
    def __init__(self, min_seed_length=None):
        self.seeds = None
        self.filter = None
        self.smems = None
        self.reads = None
        self.ref_pack = None
        self.min_seed_length = min_seed_length

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = get_mmi_parameter_set()
        index = libMA.MinimizerIndex(p_m, prefix + "genomes/slice_" + str(x) + ".mmi")
        self.index.set_mid_occ(local_max_ambiguity_fmd)
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        self.reads = reader.cpp_module.read_all()
        str_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            str_reads.append(str(nuc_seq))
        minimizer_seeds = index.seed(str_reads, self.ref_pack)
        lumper = libMA.SeedLumping(p_m)
        self.seeds = lumper.cpp_module.lump(minimizer_seeds, self.reads, self.ref_pack)
        del minimizer_seeds
        #extender = libMA.SeedExtender(p_m)
        #self.seeds = extender.cpp_module.extend(lumped_seeds, self.reads, self.ref_pack)
        if not self.min_seed_length is None:
            self.seeds = libMA.MinLength(p_m, self.min_seed_length).cpp_module.filter(self.seeds)
        self.filter = libMA.MaxExtendedToSMEM(p_m)

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        p_m = get_mmi_parameter_set()
        self.ref_pack = ref_pack
        self.reads = reads
        std_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            std_reads.append(str(nuc_seq))
        minimizer_seeds = mm_index.seed(std_reads, self.ref_pack)
        lumper = libMA.SeedLumping(p_m)
        self.seeds = lumper.cpp_module.lump(minimizer_seeds, self.reads, self.ref_pack)
        del minimizer_seeds
        #extender = libMA.SeedExtender(p_m)
        #self.seeds = extender.cpp_module.extend(lumped_seeds, self.reads, self.ref_pack)
        if not self.min_seed_length is None:
            self.seeds = libMA.MinLength(p_m, self.min_seed_length).cpp_module.filter(self.seeds)
        self.filter = libMA.MaxExtendedToSMEM(p_m)

    def run(self):
        self.smems = self.filter.cpp_module.filter(self.seeds)
        return self.smems

    def post(self):
        pass
        #for seed, read in zip(self.smems, self.reads):
        #    seed.confirm_seed_positions(read, self.ref_pack, True)

class MinimizerToMaxSpan:
    def __init__(self, min_seed_length=None):
        self.seeds = None
        self.filter = None
        self.smems = None
        self.reads = None
        self.ref_pack = None
        self.min_seed_length = min_seed_length

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = get_mmi_parameter_set()
        index = libMA.MinimizerIndex(p_m, prefix + "genomes/slice_" + str(x) + ".mmi")
        self.index.set_mid_occ(local_max_ambiguity_fmd)
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        self.reads = reader.cpp_module.read_all()
        str_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            str_reads.append(str(nuc_seq))
        minimizer_seeds = index.seed(str_reads, self.ref_pack)
        lumper = libMA.SeedLumping(p_m)
        lumped_seeds = lumper.cpp_module.lump(minimizer_seeds, self.reads, self.ref_pack)
        del minimizer_seeds
        extender = libMA.SeedExtender(p_m)
        self.seeds = extender.cpp_module.extend(lumped_seeds, self.reads, self.ref_pack)
        if not self.min_seed_length is None:
            self.seeds = libMA.MinLength(p_m, self.min_seed_length).cpp_module.filter(self.seeds)
        self.filter = libMA.MaxExtendedToMaxSpanning(p_m)

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        p_m = get_mmi_parameter_set()
        self.ref_pack = ref_pack
        self.reads = reads
        std_reads = libMA.StringVector()
        for nuc_seq in self.reads:
            std_reads.append(str(nuc_seq))
        minimizer_seeds = mm_index.seed(std_reads, self.ref_pack)
        lumper = libMA.SeedLumping(p_m)
        lumped_seeds = lumper.cpp_module.lump(minimizer_seeds, self.reads, self.ref_pack)
        del minimizer_seeds
        extender = libMA.SeedExtender(p_m)
        self.seeds = extender.cpp_module.extend(lumped_seeds, self.reads, self.ref_pack)
        if not self.min_seed_length is None:
            self.seeds = libMA.MinLength(p_m, self.min_seed_length).cpp_module.filter(self.seeds)
        self.filter = libMA.MaxExtendedToMaxSpanning(p_m)

    def run(self):
        self.smems = self.filter.cpp_module.filter(self.seeds)
        return self.smems

    def post(self):
        pass
        #for seed, read in zip(self.smems, self.reads):
        #    seed.confirm_seed_positions(read, self.ref_pack, True)


class ComputeMEMs:
    def __init__(self, ref_seq_filename, min_seed_length):
        self.reads_file_name = None
        self.ref_seq_filename = ref_seq_filename
        self.min_seed_length = min_seed_length
        self.x = 1

    def prep(self, x, prefix, paired=True):
        raise Exception("unimplemented")

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        self.x = x
        self.reads_file_name = reads_file_name

    def run(self):
        if self.x > 0.2:
            backwardMEM(self.ref_seq_filename, self.reads_file_name, self.min_seed_length)
        return []

    def post(self):
        pass

class ComputeMaxExtendedSeeds:
    def __init__(self, min_seed_length, do_smems=True, do_mems=False, extend_only=False):
        self.index = None
        self.reads = None
        self.seeder = None
        self.seeds = None
        self.ref_pack = None
        self.min_seed_length = min_seed_length
        self.extend_only = extend_only
        self.do_smems = do_smems
        self.do_mems = do_mems
        assert not (do_smems and do_mems)

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = ParameterSetManager()
        if self.do_smems:
            p_m.by_name("Seeding Technique").set(1)
        if self.do_mems:
            p_m.by_name("Seeding Technique").set(2)
        print("Seeding Technique =", p_m.by_name("Seeding Technique").get())
        p_m.by_name("Minimal Seed Length").set(self.min_seed_length)
        p_m.by_name("Maximal Ambiguity").set(local_max_ambiguity_fmd)
        self.index = FMIndex()
        self.index.load(prefix + "genomes/slice_" + str(x))
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        self.reads = reader.cpp_module.read_all()

        self.seeder = BinarySeeding(p_m)

    def prep_noise(self, ref_pack, reads, mm_index, fm_index, x, reads_file_name):
        p_m = ParameterSetManager()
        self.ref_pack = ref_pack
        if self.do_smems:
            p_m.by_name("Seeding Technique").set(1)
        if self.do_mems:
            p_m.by_name("Seeding Technique").set(2)
        print("Seeding Technique =", p_m.by_name("Seeding Technique").get())
        p_m.by_name("Minimal Seed Length").set(self.min_seed_length)
        p_m.by_name("Maximal Ambiguity").set(local_max_ambiguity_fmd)
        self.index = fm_index
        self.reads = reads
        if x > 0.2 or not self.do_mems:
            self.seeder = BinarySeeding(p_m)
        else:
            self.seeder = None

    def run(self):
        if self.seeder is None:
            return None
        if self.extend_only:
            self.seeder.cpp_module.get_segments(self.index, self.reads)
            #segsvec = self.seeder.cpp_module.get_segments(self.index, self.reads)
            #cnt = 0
            #for segs in segsvec:
            #    for seg in segs:
            #        if seg.sa_size() > 100:
            #            print(seg.sa_size(), end=" ")
            #        else:
            #            cnt += 1
            #print()
            #print("ommited:", cnt)
        else:
            self.seeds = self.seeder.cpp_module.seed(self.index, self.reads)
        return self.seeds

    def post(self):
        x = 0
        if not self.seeds is None:
            for seeds in self.seeds:
                x += len(seeds)
        print("num seeds:", x)
        pass
        #for seed, read in zip(self.seeds, self.reads):
        #    seed.confirm_seed_positions(read, self.ref_pack, True)

class CreateFmdIndex:
    def __init__(self):
        self.genome_prefix = ""
        self.pack = None
        self.fm_index = None

    def prep(self, x, prefix, paired=True):
        self.genome_prefix = prefix + "genomes/slice_" + str(x)
        self.pack = Pack()
        self.pack.load(self.genome_prefix)

    def run(self):
        self.fm_index = FMIndex(self.pack)

    def post(self):
        self.fm_index.store(self.genome_prefix)


def compareSeedSets(prefix, log_file_name, min_seed_length, do_smems=True,
                    with_single=True, with_paired=False, time_steps=linear_range):
    if x_axis_unit == "read_noise":
        ref_pack = Pack()
        ref_pack.load(reference_genome_path)
        p_m = get_mmi_parameter_set()
        mm_index = libMA.MinimizerIndex(p_m, prefix + "genomes/full.mmi")
        mm_index.set_mid_occ(local_max_ambiguity_fmd)
        fm_index = FMIndex()
        fm_index.load(reference_genome_path)
    if do_smems:
        mmi_caller = MinimizerToSmem()
    else:
        mmi_caller = MinimizerToMaxSpan()
    smem_caller = ComputeMaxExtendedSeeds(min_seed_length=min_seed_length, do_smems=do_smems)
    print(log_file_name, ": computing difference between seed sets...")
    if with_paired:
        with open(prefix + "illumina_" + log_file_name, "w") as log_file:
            log_file.write("index\tgenome size\tunique mmis\tunique smems\tshared seeds\t" +
                        "unique mmi entropy\tunique smems entropy\tshared seeds entropy\n")
            print("unique_mmis", "shared_seeds", "unique_smems", "%unique_mmis", "%unique_smems", "error_rate")
            for idx, x in enumerate(time_steps()):
                print(idx, x, "...")
                if x_axis_unit == "genome_section_size":
                    mmi_caller.prep(x, prefix)
                    smem_caller.prep(x, prefix)
                if x_axis_unit == "read_noise":
                    reads_file_name = prefix + "reads/noise_" + str(x) + ".bwa.read1.fastq.gz"
                    reads = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads_file_name)])) \
                                .cpp_module.read_all()
                    mmi_caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                    smem_caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                mmi_seeds_vec = mmi_caller.run()
                smem_seeds_vec = smem_caller.run()
                unique_mmis = libMA.seedVector()
                shared_seeds = libMA.seedVector()
                unique_smems = libMA.seedVector()
                unique_mmis_cnt = 0
                shared_seeds_cnt = 0
                unique_smems_cnt = 0
                for mmi_seeds, smem_seeds in zip(mmi_seeds_vec, smem_seeds_vec):
                    unique_mmis_, shared_seeds_, unique_smems_ = mmi_seeds.split_seed_sets(smem_seeds)
                    unique_mmis.append(unique_mmis_)
                    unique_smems.append(unique_smems_)
                    shared_seeds.append(shared_seeds_)
                    unique_mmis_cnt += len(unique_mmis_)
                    unique_smems_cnt += len(unique_smems_)
                    shared_seeds_cnt += len(shared_seeds_)

                print(unique_mmis_cnt, shared_seeds_cnt, unique_smems_cnt,
                    100*unique_mmis_cnt/max(shared_seeds_cnt + unique_mmis_cnt, 1),
                    100*unique_smems_cnt/max(shared_seeds_cnt + unique_smems_cnt, 1),
                    x)
                log_file.write(str(idx) + "\t" + str(x) + "\t" + str(unique_mmis_cnt) + "\t" + str(unique_smems_cnt) +
                                "\t" + str(shared_seeds_cnt) +
                                "\t" + str( get_seed_entropy(unique_mmis, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\t" + str( get_seed_entropy(unique_smems, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\t" + str( get_seed_entropy(shared_seeds, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\n")
    if with_single:
        with open(prefix + "pacb_" + log_file_name, "w") as log_file:
            log_file.write("index\tgenome size\tunique mmis\tunique smems\tshared seeds\t" +
                        "unique mmi entropy\tunique smems entropy\tshared seeds entropy\n")
            print("unique_mmis", "shared_seeds", "unique_smems", "%unique_mmis", "%unique_smems", "error_rate")
            for idx, x in enumerate(time_steps()):
                print(idx, x, "...")
                if x_axis_unit == "genome_section_size":
                    mmi_caller.prep(x, prefix, paired=False)
                    smem_caller.prep(x, prefix, paired=False)
                if x_axis_unit == "read_noise":
                    reads_file_name = prefix + "reads/noise_" + str(x) + ".fasta"
                    reads = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads_file_name)])) \
                                .cpp_module.read_all()
                    mmi_caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                    smem_caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                mmi_seeds_vec = mmi_caller.run()
                smem_seeds_vec = smem_caller.run()
                unique_mmis = libMA.seedVector()
                shared_seeds = libMA.seedVector()
                unique_smems = libMA.seedVector()
                unique_mmis_cnt = 0
                shared_seeds_cnt = 0
                unique_smems_cnt = 0
                for mmi_seeds, smem_seeds, in zip(mmi_seeds_vec, smem_seeds_vec):
                    unique_mmis_, shared_seeds_, unique_smems_ = mmi_seeds.split_seed_sets(smem_seeds)
                    unique_mmis.append(unique_mmis_)
                    unique_smems.append(unique_smems_)
                    shared_seeds.append(shared_seeds_)
                    unique_mmis_cnt += len(unique_mmis_)
                    unique_smems_cnt += len(unique_smems_)
                    shared_seeds_cnt += len(shared_seeds_)

                print(unique_mmis_cnt, shared_seeds_cnt, unique_smems_cnt,
                    100*unique_mmis_cnt/max(shared_seeds_cnt + unique_mmis_cnt, 1),
                    100*unique_smems_cnt/max(shared_seeds_cnt + unique_smems_cnt, 1),
                    x)
                log_file.write(str(idx) + "\t" + str(x) + "\t" + str(unique_mmis_cnt) + "\t" + str(unique_smems_cnt) +
                                "\t" + str(shared_seeds_cnt) +
                                "\t" + str( get_seed_entropy(unique_mmis, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\t" + str( get_seed_entropy(unique_smems, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\t" + str( get_seed_entropy(shared_seeds, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\n")

    print("done")

def seed_set_entropy(caller, prefix, log_file_name, with_paired=True, with_single=True, time_steps=log_range):
    if x_axis_unit == "read_noise":
        ref_pack = Pack()
        ref_pack.load(reference_genome_path)
        p_m = get_mmi_parameter_set()
        mm_index = libMA.MinimizerIndex(p_m, prefix + "genomes/full.mmi")
        mm_index.set_mid_occ(local_max_ambiguity_fmd)
        fm_index = FMIndex()
        fm_index.load(reference_genome_path)
    print(log_file_name, ": computing entropy of seed set...")
    if with_paired:
        with open(prefix + "illumina_" + log_file_name, "w") as log_file:
            log_file.write("index\tgenome size\tentropy\tpercent read covered on ref\tratio\tpercent seeds in area\n")
            for idx, x in enumerate(time_steps()):
                print(idx, "...")
                if x_axis_unit == "genome_section_size":
                    caller.prep(x, prefix)
                if x_axis_unit == "read_noise":
                    reads_file_name = prefix + "reads/noise_" + str(x) + ".bwa.read1.fastq.gz"
                    reads = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads_file_name)])) \
                                .cpp_module.read_all()
                    caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                
                start = time.time()
                seeds = caller.run()
                end = time.time()
                #q, n = get_seed_entropy(seeds, caller.reads, caller.ref_pack, True)
                n, q = get_seed_entropy(seeds, caller.reads, caller.ref_pack, False, return_seeds_in_area=True)
                print("entropy:", q)
                log_file.write(str(idx) + "\t" + str(x) + "\t" + str(q) + "\t" + str(0) + "\t" +
                               str(end - start) + "\t" + str(n) + "\n")
    if with_single:
        with open(prefix + "pacb_" + log_file_name, "w") as log_file:
            log_file.write("index\tgenome size\tentropy\tpercent read covered on ref\tratio\tpercent seeds in area\n")
            for idx, x in enumerate(time_steps()):
                print(idx, "...")
                if x_axis_unit == "genome_section_size":
                    caller.prep(x, prefix, paired=False)
                if x_axis_unit == "read_noise":
                    reads_file_name = prefix + "reads/noise_" + str(x) + ".fasta"
                    reads = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads_file_name)])) \
                                .cpp_module.read_all()
                    caller.prep_noise(ref_pack, reads, mm_index, fm_index, x, reads_file_name)
                    
                start = time.time()
                seeds = caller.run()
                end = time.time()
                #q, n = get_seed_entropy(seeds, caller.reads, caller.ref_pack, True)
                n, q = get_seed_entropy(seeds, caller.reads, caller.ref_pack, False, name=log_file_name,
                                        return_seeds_in_area=True)
                print("entropy:", q)
                log_file.write(str(idx) + "\t" + str(x) + "\t" + str(q) + "\t" + str(0) + "\t" +
                               str(end - start) + "\t" + str(n) + "\n")

    print("done")

def render_seed_set_comp(title, prefix, file_name, names,
                         divide_y_by=1, y_range=None):
    xAxes = []
    yAxes = []
    yAxesRecall = []
    yAxesPrecision = []
    legends = []
    colors = []
    xAxes2 = []
    yAxes2 = []
    legends2 = []
    colors2 = []
    with open(prefix + file_name, "r") as tsv_file:
        lines = [x[:-1].split("\t") for x in tsv_file.readlines()]
        header = {y:x for x,y in enumerate(lines[0])}
        xAxes.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes.append([ int(x[header["unique mmis"]])/divide_y_by for x in lines[1:] ])
        yAxesRecall.append([ float(x[header["shared seeds"]]) /
                ( float(x[header["unique smems"]]) + float(x[header["shared seeds"]]) ) for x in lines[1:] ])
        yAxesPrecision.append([ float(x[header["shared seeds"]]) /
                ( float(x[header["unique mmis"]]) + float(x[header["shared seeds"]]) ) for x in lines[1:] ])
        legends.append(names[0])
        colors.append(color_scheme("red"))

        xAxes2.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes2.append([ float(x[header["unique mmi entropy"]]) for x in lines[1:] ])
        legends2.append(names[0])
        colors2.append(color_scheme("red"))

        xAxes.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes.append([ int(x[header["unique smems"]])/divide_y_by for x in lines[1:] ])
        legends.append(names[1])
        colors.append(color_scheme("blue"))

        xAxes2.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes2.append([ float(x[header["unique smems entropy"]]) for x in lines[1:] ])
        legends2.append(names[1])
        colors2.append(color_scheme("blue"))
        
        xAxes.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes.append([ int(x[header["shared seeds"]])/divide_y_by for x in lines[1:] ])
        legends.append(names[2])
        colors.append(color_scheme("green"))

        xAxes2.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes2.append([ float(x[header["shared seeds entropy"]]) for x in lines[1:] ])
        legends2.append(names[2])
        colors2.append(color_scheme("green"))

    reset_output()
    if x_axis_unit == "genome_section_size":
        plot = figure(title=title, x_axis_type="log", y_range=(0,1))
    if x_axis_unit == "read_noise":
        plot = figure(title=title, y_range=(0,1))
    for yAxis, xAxis in zip(yAxesRecall, xAxes):
        plot.line(x=xAxis, y=yAxis, legend_label="recall " + file_name, line_color=color_scheme("blue"), line_width=point_to_px(4))
        plot.x(x=xAxis, y=yAxis, legend_label="recall " + file_name, line_color=color_scheme("blue"),
                size=point_to_px(7), line_width=point_to_px(2))
    for yAxis, xAxis in zip(yAxesPrecision, xAxes):
        plot.line(x=xAxis, y=yAxis, legend_label="precision " + file_name, line_color=color_scheme("red"),
                    line_width=point_to_px(4))
        plot.x(x=xAxis, y=yAxis, legend_label="precision " + file_name, line_color=color_scheme("red"),
                size=point_to_px(7), line_width=point_to_px(2))
    plot.legend.location = "top_left"
    if x_axis_unit == "genome_section_size":
        plot.xaxis.axis_label = "section length [nt]"
    if x_axis_unit == "read_noise":
        plot.xaxis.axis_label = "read noise"
    plot.yaxis.axis_label = "recall & precision"
    style_plot(plot)
    show(plot)
    if save_plots:
        plot.output_backend = "svg"
        export_svgs(plot, filename=prefix + "svg/recall-precision-" + file_name + ".svg")

    # output recall and precesion
    reset_output()
    if x_axis_unit == "genome_section_size":
        plot = figure(title=title, x_axis_type="log", y_range=y_range)
    if x_axis_unit == "read_noise":
        plot = figure(title=title, y_range=y_range)
    for yAxis, xAxis, legend, color in zip(yAxes, xAxes, legends, colors):
        plot.line(x=xAxis, y=yAxis, legend_label=legend, line_color=color, line_width=point_to_px(4))
        plot.x(x=xAxis, y=yAxis, legend_label=legend, line_color=color, size=point_to_px(7), line_width=point_to_px(2))
        #plot.x(x=xAxis, y=yAxis, legend_label=legend, color=color)
    plot.legend.location = "top_left"
    if x_axis_unit == "genome_section_size":
        plot.xaxis.axis_label = "section length [nt]"
    if x_axis_unit == "read_noise":
        plot.xaxis.axis_label = "read noise"
    plot.yaxis.axis_label = "number of seeds"
    style_plot(plot)
    show(plot)
    if save_plots:
        plot.output_backend = "svg"
        export_svgs(plot, filename=prefix + "svg/" + file_name + ".svg")
    #output entropy as well

    reset_output()
    if x_axis_unit == "genome_section_size":
        plot = figure(title=title, x_axis_type="log")
    if x_axis_unit == "read_noise":
        plot = figure(title=title, y_range=(-1,30))
    for yAxis, xAxis, legend, color in zip(yAxes2, xAxes2, legends2, colors2):
        plot.line(x=xAxis, y=yAxis, legend_label=legend, line_color=color, line_width=point_to_px(4))
        plot.x(x=xAxis, y=yAxis, legend_label=legend, line_color=color, size=point_to_px(7), line_width=point_to_px(2))
        #plot.x(x=xAxis, y=yAxis, legend_label=legend, color=color)
    plot.legend.location = "top_left"
    if x_axis_unit == "genome_section_size":
        plot.xaxis.axis_label = "section length [nt]"
    if x_axis_unit == "read_noise":
        plot.xaxis.axis_label = "read noise"
    plot.yaxis.axis_label = "seed entropy"
    style_plot(plot)
    show(plot)
    if save_plots:
        plot.output_backend = "svg"
        export_svgs(plot, filename=prefix + "svg/entropy-" + file_name + ".svg")

def read_generation(time_steps=linear_range, backwards=False, survivor_error_profile=survivor_error_profile):
    ref_pack = Pack()
    ref_pack.load(reference_genome_path)
    print("genome size", ref_pack.unpacked_size())
    if True:
        if x_axis_unit == "genome_section_size":
            generate_genomes(time_steps, ref_pack, prefix + "genomes", backwards=backwards)
        if x_axis_unit == "read_noise":
            generate_profiles(time_steps, survivor_error_profile=survivor_error_profile)
            libMA.MinimizerIndex(get_mmi_parameter_set(), ref_pack.contigSeqs(),
                                ref_pack.contigNames()).dump(prefix + "genomes/full.mmi")
    generate_reads(time_steps, prefix + "reads/", prefix + "genomes")

def seed_entropy_analysis(w_paired=True, w_single=True, time_steps=linear_range):
    if True:
        seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small), prefix, "smem_seed_entropy.tsv",
                                                w_paired, w_single, time_steps)
        seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large), prefix,
                        "smem_"+str_msl+"_seed_entropy.tsv", w_paired, w_single, time_steps)
        seed_set_entropy(ExtendMinimizers(min_seed_length=mem_size_large), prefix, "mem_seed_entropy_"+str_msl+".tsv",
                                            w_paired, w_single, time_steps)
    if True:
        seed_set_entropy(MinimizerToSmem(min_seed_length=mem_size_small), prefix, "mmi_to_smem_seed_entropy.tsv",
                                        w_paired, w_single, time_steps)
        seed_set_entropy(MinimizerToMaxSpan(min_seed_length=mem_size_small), prefix,
                                            "mini_to_max_sp_seed_entropy.tsv", w_paired, w_single, time_steps)
    if True:
        seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False),
                        prefix, "max_sp_seed_entropy.tsv", w_paired, w_single, time_steps)
        seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large, do_smems=False), prefix,
                                               "max_sp_"+str_msl+"_seed_entropy.tsv", w_paired, w_single, time_steps)
        
        seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False, do_mems=True), prefix,
                    "fmd_mem_seed_entropy.tsv", w_paired, w_single, time_steps)
        seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large, do_smems=False, do_mems=True), prefix,
                    "fmd_mem_"+str_msl+"_seed_entropy.tsv", w_paired, w_single, time_steps)
    if True:
        seed_set_entropy(ComputeMinimizers(), prefix, "minimizer_seed_entropy.tsv", w_paired, w_single, time_steps)
        seed_set_entropy(ExtendMinimizers(), prefix, "mem_seed_entropy.tsv", w_paired, w_single, time_steps)
    if True:
        render_times("Seed entropy - illumina", prefix, [
                    #[("MEMs l≥19", "illumina_fmd_mem_seed_entropy.tsv", "grey")],
                    #[("MEMs l≥"+str_msl, "illumina_fmd_mem_"+str_msl+"_seed_entropy.tsv", "black")],
                    #[("SMEMs l≥"+str_mss, "illumina_smem_seed_entropy.tsv", "lightblue")],
                    #[("SMEMs l≥"+str_msl, "illumina_smem_"+str_msl+"_seed_entropy.tsv", "blue")],
                    #[("maximally spanning l≥"+str_mss, "illumina_max_sp_seed_entropy.tsv", "lightgreen")],
                    #[("maximally spanning l≥"+str_msl+"", "illumina_max_sp_"+str_msl+"_seed_entropy.tsv", "green")],
                    [("Alg. 2a l≥"+str_mss+" (SMEMs)", "illumina_mmi_to_smem_seed_entropy.tsv", "pink")],
                    [("Alg. 1 l≥"+str_mss+" (MEMs)", "illumina_mem_seed_entropy.tsv", "orange")],
                    #[("Alg. 1 l≥"+str_msl+" (MEMs)", "illumina_mem_seed_entropy_"+str_msl+".tsv", "yellow")],
                    [(str_mm+"-minimizer", "illumina_minimizer_seed_entropy.tsv", "red")],
                    [("Alg. 2b l≥"+str_mss+" (max. spanning)", "illumina_mini_to_max_sp_seed_entropy.tsv", "purple")],
                ], 
                "illumina_seed_entropy",
                yAxisKey="entropy",
                #yAxisKey2="percent read covered on ref",
                y_axis_log=True )
    if True:
        render_times("Seed entropy - PacBio", prefix, [
                            #[("MEMs l≥"+str_mss, "pacb_fmd_mem_seed_entropy.tsv", "grey")],
                            #[("MEMs l≥"+str_msl, "pacb_fmd_mem_"+str_msl+"_seed_entropy.tsv", "black")],
                            #[("SMEMs l≥"+str_mss, "pacb_smem_seed_entropy.tsv", "lightblue")],
                            #[("SMEMs l≥"+str_msl, "pacb_smem_"+str_msl+"_seed_entropy.tsv", "blue")],
                            #[("maximally spanning l≥"+str_mss, "pacb_max_sp_seed_entropy.tsv", "lightgreen")],
                            #[("maximally spanning l≥"+str_msl, "pacb_max_sp_"+str_msl+"_seed_entropy.tsv", "green")],
                            [("Alg. 2a l≥"+str_mss+" (SMEMs)", "pacb_mmi_to_smem_seed_entropy.tsv", "pink")],
                            [("Alg. 1 l≥"+str_mss+" (MEMs)", "pacb_mem_seed_entropy.tsv", "orange")],
                            #[("Alg. 1 l≥"+str_msl+" (MEMs)", "pacb_mem_seed_entropy_"+str_msl+".tsv", "yellow")],
                            [(str_mm+"-minimizer", "pacb_minimizer_seed_entropy.tsv", "red")],
                            [("Alg. 2b l≥"+str_mss+" (max. spanning)", "pacb_mini_to_max_sp_seed_entropy.tsv", "purple")],
                        ],
                        "pacb_seed_entropy",
                        yAxisKey="entropy",
                        #yAxisKey2="percent read covered on ref",
                        y_axis_log=True )
        render_times("Percent Seeds in correct area - PacBio", prefix, [
                            #[("MEMs l≥"+str_mss, "pacb_fmd_mem_seed_entropy.tsv", "grey")],
                            #[("MEMs l≥"+str_msl, "pacb_fmd_mem_"+str_msl+"_seed_entropy.tsv", "black")],
                            #[("SMEMs l≥"+str_mss, "pacb_smem_seed_entropy.tsv", "lightblue")],
                            #[("SMEMs l≥"+str_msl, "pacb_smem_"+str_msl+"_seed_entropy.tsv", "blue")],
                            #[("maximally spanning l≥"+str_mss, "pacb_max_sp_seed_entropy.tsv", "lightgreen")],
                            #[("maximally spanning l≥"+str_msl, "pacb_max_sp_"+str_msl+"_seed_entropy.tsv", "green")],
                            [("Alg. 2a l≥"+str_mss+" (SMEMs)", "pacb_mmi_to_smem_seed_entropy.tsv", "pink")],
                            [("Alg. 1 l≥"+str_mss+" (MEMs)", "pacb_mem_seed_entropy.tsv", "orange")],
                            #[("Alg. 1 l≥"+str_msl+" (MEMs)", "pacb_mem_seed_entropy_"+str_msl+".tsv", "yellow")],
                            [(str_mm+"-minimizer", "pacb_minimizer_seed_entropy.tsv", "red")],
                            [("Alg. 2b l≥"+str_mss+" (max. spanning)", "pacb_mini_to_max_sp_seed_entropy.tsv", "purple")],
                        ],
                        "pacb_seed_correct_rate",
                        yAxisKey="percent seeds in area",
                        #yAxisKey2="percent read covered on ref",
                        y_axis_log=True )
    if True:
        render_times("Seed entropy ratio - PacBio", prefix, [
                            #[("MEMs l≥"+str_mss, "pacb_fmd_mem_seed_entropy.tsv", "grey")],
                            #[("MEMs l≥"+str_msl, "pacb_fmd_mem_"+str_msl+"_seed_entropy.tsv", "black")],
                            #[("SMEMs l≥"+str_mss, "pacb_smem_seed_entropy.tsv", "lightblue")],
                            #[("SMEMs l≥"+str_msl, "pacb_smem_"+str_msl+"_seed_entropy.tsv", "blue")],
                            #[("maximally spanning l≥"+str_mss, "pacb_max_sp_seed_entropy.tsv", "lightgreen")],
                            #[("maximally spanning l≥"+str_msl, "pacb_max_sp_"+str_msl+"_seed_entropy.tsv", "green")],
                            [("Alg. 2a l≥"+str_mss+" (SMEMs)", "pacb_mmi_to_smem_seed_entropy.tsv", "pink")],
                            [("Alg. 1 l≥"+str_mss+" (MEMs)", "pacb_mem_seed_entropy.tsv", "orange")],
                            #[("Alg. 1 l≥"+str_msl+" (MEMs)", "pacb_mem_seed_entropy_"+str_msl+".tsv", "yellow")],
                            [(str_mm+"-minimizer", "pacb_minimizer_seed_entropy.tsv", "red")],
                            [("Alg. 2b l≥"+str_mss+" (max. spanning)", "pacb_mini_to_max_sp_seed_entropy.tsv", "purple")],
                        ],
                        "pacb_seed_entropy_ratio",
                        yAxisKey="entropy",
                        yAxisKeyDivident="ratio",
                        #yAxisKey2="percent read covered on ref",
                        y_axis_log=True )



def runtime_analysis():
    if x_axis_unit == "genome_section_size":
        measure_time(CreateFmdIndex(), prefix, "fm_index_construction.tsv")
        measure_time(CreateMinimizerIndex(), prefix, "minimizer_index_construction.tsv")
    if True:
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small), prefix, "smem_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, extend_only=True), prefix,
                                            "smem_extend_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False), prefix,
                    "max_sp_seed_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False, extend_only=True), prefix,
                    "max_sp_seed_extend_computation.tsv", True)
    if False: # crashes on full genome...
        measure_time(ComputeMEMs(reference_genome_fasta, mem_size_large), prefix,
                    "mem_"+str_msl+"_seed_computation.tsv", True)
        measure_time(ComputeMEMs(reference_genome_fasta, mem_size_small), prefix,
                    "mem_seed_computation.tsv", True)
    if True:
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False, do_mems=True), prefix,
                    "fmd_mem_seed_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False, do_mems=True,
                                             extend_only=True), prefix,
                    "fmd_mem_seed_extend_computation.tsv", True)

        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large), prefix,
                    "smem_"+str_msl+"_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large, extend_only=True), prefix,
                    "smem_extend_"+str_msl+"_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(do_smems=False, min_seed_length=mem_size_large), prefix,
                    "max_sp_"+str_msl+"_seed_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(do_smems=False, min_seed_length=mem_size_large, extend_only=True), prefix,
                    "max_sp_"+str_msl+"_seed_extend_computation.tsv", True)
    if True:
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large, do_smems=False, do_mems=True), prefix,
                    "fmd_mem_"+str_msl+"_seed_computation.tsv", True)
        measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large, do_smems=False, do_mems=True, 
                                             extend_only=True), prefix,
                    "fmd_mem_"+str_msl+"_seed_extend_computation.tsv", True)
    if True:
        measure_time(ComputeMinimizers(), prefix, "minimizer_computation.tsv", True)
        measure_time(LumpMinimizers(), prefix, "minimizer_lumping.tsv", True)
    if True:
        measure_time(MinimizerToSmem(), prefix, "minimizer_to_smem.tsv", True)
        measure_time(MinimizerToMaxSpan(), prefix, "minimizer_to_max_span.tsv", True)
        measure_time(ExtendThenSortMinimizers(), prefix, "extend_then_sort.tsv", True)

    if False:
        render_times("runtimes - illumina", prefix, [
                        #[("MEMs l≥"+str_mss, "illumina_mem_seed_computation.tsv", "black")],
                        #[("MEMs l≥"+str_msl, "illumina_mem_"+str_msl+"_seed_computation.tsv", "grey")],
                        #[("SMEMs l≥"+str_mss, "illumina_smem_computation.tsv", "blue")],
                        #[("SMEMs l≥"+str_msl, "illumina_smem_"+str_msl+"_computation.tsv", "lightblue")],
                        #[("maximally spanning l≥"+str_mss, "illumina_max_sp_seed_computation.tsv", "green")],
                        #[("maximally spanning l≥"+str_msl, "illumina_max_sp_"+str_msl+"_seed_computation.tsv",
                        #"lightgreen")],
        
                        #[(str_mm+"-minimizer", "illumina_minimizer_computation.tsv", "red"),
                        #("Alg. 1", "illumina_minimizer_lumping.tsv", "orange"),
                        #("Alg. 2a (SMEMs)", "illumina_minimizer_to_smem.tsv", "pink")],
        
                        [#(str_mm+"-minimizer", "illumina_minimizer_computation.tsv", "red"),
                        ("Alg. 1", "illumina_minimizer_lumping.tsv", "orange")
                        #,("Alg. 2b (max. spanning)", "illumina_minimizer_to_max_span.tsv", "purple")
                        ],
        
                        [#(str_mm+"-minimizer", "illumina_minimizer_computation.tsv", "red"),
                        ("Extend-Filter", "illumina_extend_then_sort.tsv", "yellow")],
                    ],
                    "illumina_times",
                    divide_y_by=num_illumina_reads/1000 )
    if True:
        render_times("runtimes - illumina", prefix, [
                        [("MEMs l≥"+str_mss, "illumina_fmd_mem_seed_computation.tsv", "black")],
                        [("MEMs l≥"+str_msl, "illumina_fmd_mem_"+str_msl+"_seed_computation.tsv", "grey")],
                        [("SMEMs l≥"+str_mss, "illumina_smem_computation.tsv", "blue")],
                        [("SMEMs l≥"+str_msl, "illumina_smem_"+str_msl+"_computation.tsv", "lightblue")],
                        [("maximally spanning l≥"+str_mss, "illumina_max_sp_seed_computation.tsv", "green")],
                        [("maximally spanning l≥"+str_msl, "illumina_max_sp_"+str_msl+"_seed_computation.tsv",
                        "lightgreen")],
        
                        [(str_mm+"-minimizer", "illumina_minimizer_computation.tsv", "red"),
                        ("Alg. 1", "illumina_minimizer_lumping.tsv", "orange"),
                        ("Alg. 2a (SMEMs)", "illumina_minimizer_to_smem.tsv", "pink")],
        
                        [(str_mm+"-minimizer", "illumina_minimizer_computation.tsv", "red"),
                        ("Alg. 1", "illumina_minimizer_lumping.tsv", "orange"),
                        ("Alg. 2b (max. spanning)", "illumina_minimizer_to_max_span.tsv", "purple")],
        
                        [(str_mm+"-minimizer", "illumina_minimizer_computation.tsv", "red"),
                        ("Extend-Filter", "illumina_extend_then_sort.tsv", "yellow")],
                    ],
                    "illumina_times_2",
                    divide_y_by=num_illumina_reads/1000, y_axis_log=True, y_range=(0.02,11) )
    l = [
        #[("fm-MEMs l≥"+str_mss, "pacb_mem_seed_computation.tsv", "black")],
        #[("fm-MEMs l≥"+str_msl, "pacb_mem_"+str_msl+"_seed_computation.tsv", "grey")],
        [("MEMs l≥"+str_mss, "pacb_fmd_mem_seed_computation.tsv", "black")],
        [("MEMs l≥"+str_msl, "pacb_fmd_mem_"+str_msl+"_seed_computation.tsv", "grey")],
        [("MEMs l≥"+str_mss, "pacb_fmd_mem_seed_extend_computation.tsv", "black")],
        [("MEMs l≥"+str_msl, "pacb_fmd_mem_"+str_msl+"_seed_extend_computation.tsv", "grey")],

        [("SMEMs l≥"+str_mss, "pacb_smem_computation.tsv", "blue")],
        [("SMEMs l≥"+str_msl, "pacb_smem_"+str_msl+"_computation.tsv", "lightblue")],
        [("SMEMs l≥"+str_mss, "pacb_smem_extend_computation.tsv", "blue")],
        [("SMEMs l≥"+str_msl, "pacb_smem_extend_"+str_msl+"_computation.tsv", "lightblue")],

        [("maximally spanning l≥"+str_mss, "pacb_max_sp_seed_computation.tsv", "green")],
        [("maximally spanning l≥"+str_msl, "pacb_max_sp_"+str_msl+"_seed_computation.tsv", "lightgreen")],
        [("maximally spanning l≥"+str_mss, "pacb_max_sp_seed_extend_computation.tsv", "green")],
        [("maximally spanning l≥"+str_msl, "pacb_max_sp_"+str_msl+"_seed_extend_computation.tsv", "lightgreen")],

        [(str_mm+"-minimizer", "pacb_minimizer_computation.tsv", "red"),
        ("Alg. 1", "pacb_minimizer_lumping.tsv", "orange")
        ,("Alg. 2a (SMEMs)", "pacb_minimizer_to_smem.tsv", "pink")
        ],

        [(str_mm+"-minimizer", "pacb_minimizer_computation.tsv", "red"),
        ("Alg. 1", "pacb_minimizer_lumping.tsv", "orange"),
        ("Alg. 2b (max. spanning)", "pacb_minimizer_to_max_span.tsv", "purple")],

        [(str_mm+"-minimizer", "pacb_minimizer_computation.tsv", "red"),
        ("Extend-Filter", "pacb_extend_then_sort.tsv", "yellow")],
    ]
    if True:
        render_times("runtimes - PacBio", prefix, l, "pacb_times",
                    divide_y_by=num_pacb_reads/1000, y_range=(0.5,600), y_axis_log=True )
    if False:
        render_times("runtimes - PacBio - fmd vs lcp", prefix, [
                        [("LCP-MEMs l≥"+str_mss, "pacb_mem_seed_computation.tsv", "blue")],
                        [("LCP-MEMs l≥"+str_msl, "pacb_mem_"+str_msl+"_seed_computation.tsv", "lightblue")],
                        [("FMD-MEMs l≥"+str_mss, "pacb_fmd_mem_seed_computation.tsv", "green")],
                        [("FMD-MEMs l≥"+str_msl, "pacb_fmd_mem_"+str_msl+"_seed_computation.tsv", "lightgreen")],
                    ], "pacb_times_fm_lcp",
                    divide_y_by=num_pacb_reads/1000, y_axis_log=True,
                    y_range=(70,2800) )
    #render_times("runtime - index creation", prefix, [
    #                [("FMD-Index", "illumina_fm_index_construction.tsv", "blue")],
    #                [("10,19-minimizer Index", "illumina_minimizer_index_construction.tsv", "red")],
    #            ],
    #            "index_times",
    #            y_axis_log=True,
    #            divide_y_by=0.001 )

def seed_set_diff_analysis():
    if True:
        compareSeedSets(prefix, "seed_set_comparison.tsv", mem_size_small, with_single=True)
        compareSeedSets(prefix, "seed_set_comparison_"+str_msl+".tsv", mem_size_large, with_single=True)
        compareSeedSets(prefix, "seed_set_comparison_max_sp.tsv", mem_size_small, do_smems=False, with_single=True)
        compareSeedSets(prefix, "seed_set_comparison_max_sp_"+str_msl+".tsv", mem_size_large, do_smems=False,
                        with_single=True)

    if False:
        # we do not actually compute this analysis for illumina reads - so we cannot render it
        render_seed_set_comp("seed set comparison - illumina", prefix, "illumina_seed_set_comparison.tsv",
                            ("unique "+str_mm+"-minimizers", "unique SMEMs l≥"+str_mss, "shared seeds l≥"+str_mss),
                            divide_y_by=num_illumina_reads, y_range=(0,11))
        render_seed_set_comp("seed set comparison - illumina", prefix, "illumina_seed_set_comparison_"+str_msl+".tsv",
                            ("unique "+str_mm+"-minimizers", "unique SMEMs l≥"+str_msl+"", "shared seeds l≥"+str_msl), divide_y_by=num_illumina_reads, y_range=(0,11))
        render_seed_set_comp("seed set comparison - illumina", prefix,
                            "illumina_seed_set_comparison_max_sp.tsv",
                            ("unique "+str_mm+"-minimizers", "unique max. spanning l≥"+str_mss, "shared seeds l≥"+str_mss),
                            divide_y_by=num_illumina_reads, y_range=(0,11))
        render_seed_set_comp("seed set comparison - illumina", prefix,
                            "illumina_seed_set_comparison_max_sp_"+str_msl+".tsv",
                            ("unique "+str_mm+"-minimizers", "unique max. spanning l≥"+str_msl, "shared seeds l≥"+str_msl),
                            divide_y_by=num_illumina_reads, y_range=(0,11))

    if True:
        render_seed_set_comp("seed set comparison - pacBio", prefix, "pacb_seed_set_comparison.tsv",
                            ("unique "+str_mm+"-minimizers", "unique SMEMs l≥"+str_mss, "shared seeds l≥"+str_mss),
                            divide_y_by=num_pacb_reads, y_range=(-42, 922))
        render_seed_set_comp("seed set comparison - pacBio", prefix, "pacb_seed_set_comparison_"+str_msl+".tsv",
                            ("unique "+str_mm+"-minimizers", "unique SMEMs l≥"+str_msl, "shared seeds l≥"+str_msl), divide_y_by=num_pacb_reads, y_range=(-42, 922))
    if True:
        render_seed_set_comp("seed set comparison - pacBio", prefix, "pacb_seed_set_comparison_max_sp.tsv",
                            ("unique "+str_mm+"-minimizers", "unique max. spanning l≥"+str_mss, "shared seeds l≥"+str_mss),
                            divide_y_by=num_pacb_reads, y_range=(0, 800))
        render_seed_set_comp("seed set comparison - pacBio", prefix,
                            "pacb_seed_set_comparison_max_sp_"+str_msl+".tsv",
                            ("unique "+str_mm+"-minimizers", "unique max. spanning l≥"+str_msl, "shared seeds l≥"+str_msl),
                            divide_y_by=num_pacb_reads, y_range=(0, 800))


read_generation()
runtime_analysis()
seed_entropy_analysis()
seed_set_diff_analysis()
