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

def generate_genomes(time_steps, ref_pack, out_prefix, backwards=False):
    if not os.path.exists(prefix + "svg/"):
        os.mkdir(prefix + "svg/")
    if not os.path.exists(prefix + "genomes/"):
        os.mkdir(prefix + "genomes/")
    if not os.path.exists(prefix + "reads/"):
        os.mkdir(prefix + "reads/")
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

def generate_reads(time_steps, out_prefix, genome_prefix):
    print("generating reads...")
    for idx, x in enumerate(time_steps()):
        print(idx, "...")
        dwgsim.create_illumina_reads_dwgsim(genome_prefix + "/slice_" + str(x) + ".fasta", out_prefix,
                                            num_illumina_reads, "slice_" + str(x), illumina_read_size)
        dwgsim.create_reads_survivor(genome_prefix + "/slice_" + str(x) + ".fasta", out_prefix, num_pacb_reads,
                                     "slice_" + str(x), survivor_error_profile)
    print("done")

def measure_time(caller, prefix, log_file_name, with_single=False, time_steps=log_range):
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

def render_times(title, prefix, element_list, out_file, yAxisKey="runtime", y_axis_log=False, yAxisKey2=None,
                 divide_y_by=1):
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
                xAxis_ = [ int(x[header["genome size"]]) for x in lines[1:] ]
                #assert xAxis is None or xAxis == xAxis_
                yAxes_ = [ float(x[header[yAxisKey]])/divide_y_by + y for x,y in zip(lines[1:], currY) ]
                xAxes.append(xAxis_)
                yAxes.append(yAxes_)
                currY = yAxes_
                if not yAxisKey2 is None:
                    yAxes_2 = [ float(x[header[yAxisKey2]])/divide_y_by for x in lines[1:] ]
                    yAxes.append(yAxes_2)
    reset_output()
    if y_axis_log:
        plot = figure(title=title, x_axis_type="log", y_axis_type="log")
    else:
        plot = figure(title=title, x_axis_type="log")
    step = 1 if yAxisKey2 is None else 2
    for yAxis, xAxis, legend, color in zip(yAxes[::step], xAxes, legends, colors):
        plot.line(x=xAxis, y=yAxis, legend=legend, line_color=color, line_width=point_to_px(4))
        #plot.x(x=xAxis, y=yAxis, legend=legend, color=color)
    if not yAxisKey2 is None:
        for yAxis, xAxis, color in zip(yAxes[1::2], xAxes, colors):
            plot.line(x=xAxis, y=yAxis, line_color=color, line_dash=[1,1], line_width=point_to_px(4))
            #plot.cross(x=xAxis, y=yAxis, color=color)
    plot.legend.location = "top_left"
    plot.xaxis.axis_label = "section length [nt]"
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

def get_seed_entropy(seeds_list, reads, pack, also_return_percent_covered=False):
    read_pos = get_read_positions(reads, pack)
    # accumulated overlap with read region
    hits = 0
    # accumulated overlap with non read region
    num_nuc = 0
    # accumulated read size
    r_sum = 0

    num_seeds = 0
    for read in reads:
        r_sum += len(read)
    for (read_start, read_end, is_forward_strand), seeds in zip(read_pos, seeds_list):
        seeds.sort_by_ref_pos()
        max_r_pos = 0
        num_seeds += len(seeds)
        for seed in seeds:
            seed_start = max(seed.start_ref, max_r_pos)
            seed_end = seed.start_ref + seed.size
            max_r_pos = max(max_r_pos, seed_end)

            num_nuc += overlap(0, read_start, seed_start, seed_end)
            num_nuc += overlap(read_end, pack.unpacked_size(), seed_start, seed_end)
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
            return 0
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

    def prep(self, x, prefix, paired=True):
        ref_pack = Pack()
        ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = get_mmi_parameter_set()
        index = libMA.MinimizerIndex(p_m, prefix + "genomes/slice_" + str(x) + ".mmi")
        if paired:
            reads1 = prefix + "reads/slice_" + str(x) + ".bwa.read1.fastq.gz"
            reads2 = prefix + "reads/slice_" + str(x) + ".bwa.read2.fastq.gz"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        else:
            reads1 = prefix + "reads/slice_" + str(x) + ".fasta"
            reader = libMA.FileListReader(p_m, libMA.filePathVector([libMA.path(reads1)]))
        nuc_seq_reads = reader.cpp_module.read_all()
        reads = libMA.StringVector()
        for nuc_seq in nuc_seq_reads:
            reads.append(str(nuc_seq))
        self.seeds = index.seed(reads, ref_pack)
        self.lumper = libMA.SeedLumping(p_m)

    def run(self):
        self.lumper.cpp_module.lump(self.seeds)

    def post(self):
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
        self.seeds = lumper.cpp_module.lump(minimizer_seeds)
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
        lumped_seeds = lumper.cpp_module.lump(minimizer_seeds)
        del minimizer_seeds
        extender = libMA.SeedExtender(p_m)
        self.seeds = extender.cpp_module.extend(lumped_seeds, self.reads, self.ref_pack)
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
        lumped_seeds = lumper.cpp_module.lump(minimizer_seeds)
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


class ComputeMaxExtendedSeeds:
    def __init__(self, min_seed_length, do_smems=True):
        self.index = None
        self.reads = None
        self.seeder = None
        self.seeds = None
        self.ref_pack = None
        self.min_seed_length = min_seed_length-1
        self.do_smems = do_smems

    def prep(self, x, prefix, paired=True):
        self.ref_pack = Pack()
        self.ref_pack.load(prefix + "genomes/slice_" + str(x))
        p_m = ParameterSetManager()
        if self.do_smems:
            p_m.by_name("Seeding Technique").set(1)
        print("Seeding Technique=", p_m.by_name("Seeding Technique").get())
        p_m.by_name("Minimal Seed Length").set(self.min_seed_length)
        p_m.by_name("Maximal Ambiguity").set(max_ambiguity_fmd)
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

    def run(self):
        self.seeds = self.seeder.cpp_module.seed(self.index, self.reads)
        return self.seeds

    def post(self):
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
                    with_single=False, time_steps=log_range):
    #mmi_caller = ExtendMinimizers()
    if do_smems:
        mmi_caller = MinimizerToSmem()
    else:
        mmi_caller = MinimizerToMaxSpan()
    smem_caller = ComputeMaxExtendedSeeds(min_seed_length=min_seed_length, do_smems=do_smems)
    print(log_file_name, ": computing difference between seed sets...")
    with open(prefix + "illumina_" + log_file_name, "w") as log_file:
        log_file.write("index\tgenome size\tunique mmis\tunique smems\tshared seeds\t" +
                       "unique mmi entropy\tunique smems entropy\tshared seeds entropy\n")
        print("unique_mmis", "shared_seeds", "unique_smems", "%unique_mmis", "%unique_smems")
        for idx, x in enumerate(time_steps()):
            print(idx, "...")
            mmi_caller.prep(x, prefix)
            smem_caller.prep(x, prefix)
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
                  100*unique_smems_cnt/max(shared_seeds_cnt + unique_smems_cnt, 1))
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
            print("unique_mmis", "shared_seeds", "unique_smems")
            for idx, x in enumerate(time_steps()):
                print(idx, "...")
                mmi_caller.prep(x, prefix, paired=False)
                smem_caller.prep(x, prefix, paired=False)
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
                    100*unique_smems_cnt/max(shared_seeds_cnt + unique_smems_cnt, 1))
                log_file.write(str(idx) + "\t" + str(x) + "\t" + str(unique_mmis_cnt) + "\t" + str(unique_smems_cnt) +
                                "\t" + str(shared_seeds_cnt) +
                                "\t" + str( get_seed_entropy(unique_mmis, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\t" + str( get_seed_entropy(unique_smems, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\t" + str( get_seed_entropy(shared_seeds, smem_caller.reads, smem_caller.ref_pack) ) +
                                "\n")

    print("done")

def seed_set_entropy(caller, prefix, log_file_name, with_single=False, time_steps=log_range):
    print(log_file_name, ": computing entropy of seed set...")
    with open(prefix + "illumina_" + log_file_name, "w") as log_file:
        log_file.write("index\tgenome size\tentropy\tpercent read covered on ref\n")
        for idx, x in enumerate(time_steps()):
            print(idx, "...")
            caller.prep(x, prefix)
            seeds = caller.run()
            #q, n = get_seed_entropy(seeds, caller.reads, caller.ref_pack, True)
            q = get_seed_entropy(seeds, caller.reads, caller.ref_pack, False)
            print("entropy:", q)
            log_file.write(str(idx) + "\t" + str(x) + "\t" + str(q) + "\t" + str(0) + "\n")
    if with_single:
        with open(prefix + "pacb_" + log_file_name, "w") as log_file:
            log_file.write("index\tgenome size\tentropy\tpercent read covered on ref\n")
            for idx, x in enumerate(time_steps()):
                print(idx, "...")
                caller.prep(x, prefix, paired=False)
                seeds = caller.run()
                #q, n = get_seed_entropy(seeds, caller.reads, caller.ref_pack, True)
                q = get_seed_entropy(seeds, caller.reads, caller.ref_pack, False)
                print("entropy:", q)
                log_file.write(str(idx) + "\t" + str(x) + "\t" + str(q) + "\t" + str(0) + "\n")

    print("done")

def render_seed_set_comp(title, prefix, file_name, names,
                         divide_y_by=1):
    xAxes = []
    yAxes = []
    legends = []
    colors = []
    xAxes2 = []
    yAxes2 = []
    legends2 = []
    colors2 = []
    with open(prefix + file_name, "r") as tsv_file:
        lines = [x[:-1].split("\t") for x in tsv_file.readlines()]
        header = {y:x for x,y in enumerate(lines[0])}
        xAxes.append([ int(x[header["genome size"]]) for x in lines[1:] ])
        yAxes.append([ int(x[header["unique mmis"]])/divide_y_by for x in lines[1:] ])
        legends.append(names[0])
        colors.append(color_scheme("red"))

        xAxes2.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes2.append([ float(x[header["unique mmi entropy"]]) for x in lines[1:] ])
        legends2.append(names[0])
        colors2.append(color_scheme("red"))

        xAxes.append([ int(x[header["genome size"]]) for x in lines[1:] ])
        yAxes.append([ int(x[header["unique smems"]])/divide_y_by for x in lines[1:] ])
        legends.append(names[1])
        colors.append(color_scheme("blue"))

        xAxes2.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes2.append([ float(x[header["unique smems entropy"]]) for x in lines[1:] ])
        legends2.append(names[1])
        colors2.append(color_scheme("blue"))
        
        xAxes.append([ int(x[header["genome size"]]) for x in lines[1:] ])
        yAxes.append([ int(x[header["shared seeds"]])/divide_y_by for x in lines[1:] ])
        legends.append(names[2])
        colors.append(color_scheme("green"))

        xAxes2.append([ float(x[header["genome size"]]) for x in lines[1:] ])
        yAxes2.append([ float(x[header["shared seeds entropy"]]) for x in lines[1:] ])
        legends2.append(names[2])
        colors2.append(color_scheme("green"))

    reset_output()
    plot = figure(title=title, x_axis_type="log")
    for yAxis, xAxis, legend, color in zip(yAxes, xAxes, legends, colors):
        plot.line(x=xAxis, y=yAxis, legend=legend, line_color=color, line_width=point_to_px(4))
        #plot.x(x=xAxis, y=yAxis, legend=legend, color=color)
    plot.legend.location = "top_left"
    plot.xaxis.axis_label = "section length [nt]"
    plot.yaxis.axis_label = "number of seeds"
    style_plot(plot)
    show(plot)
    if save_plots:
        plot.output_backend = "svg"
        export_svgs(plot, filename=prefix + "svg/" + file_name + ".svg")
    #output entropy as well
    reset_output()
    plot = figure(title=title, x_axis_type="log")
    for yAxis, xAxis, legend, color in zip(yAxes2, xAxes2, legends2, colors2):
        plot.line(x=xAxis, y=yAxis, legend=legend, line_color=color, line_width=point_to_px(4))
        #plot.x(x=xAxis, y=yAxis, legend=legend, color=color)
    plot.legend.location = "top_left"
    plot.xaxis.axis_label = "section length [nt]"
    plot.yaxis.axis_label = "seed entropy"
    style_plot(plot)
    show(plot)
    if save_plots:
        plot.output_backend = "svg"
        export_svgs(plot, filename=prefix + "svg/entropy-" + file_name + ".svg")

def read_generation(backwards=False):
    ref_pack = Pack()
    ref_pack.load(reference_genome_path)
    print("genome size", ref_pack.unpacked_size())
    generate_genomes(log_range, ref_pack, prefix + "genomes", backwards=backwards)
    generate_reads(log_range, prefix + "reads/", prefix + "genomes")

def seed_entropy_analysis():
    seed_set_entropy(ComputeMinimizers(), prefix, "minimizer_seed_entropy.tsv", True)
    seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small), prefix, "smem_seed_entropy.tsv", True)
    seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large), prefix,
                     "smem_"+str_msl+"_seed_entropy.tsv", True)
    seed_set_entropy(MinimizerToSmem(min_seed_length=mem_size_large), prefix, "mmi_to_smem_seed_entropy.tsv", True)
    seed_set_entropy(ExtendMinimizers(), prefix, "mem_seed_entropy.tsv", True)
    seed_set_entropy(ExtendMinimizers(min_seed_length=mem_size_large), prefix, "mem_seed_entropy_"+str_msl+".tsv", True)
    seed_set_entropy(MinimizerToMaxSpan(min_seed_length=mem_size_large), prefix,
                                        "mini_to_max_sp_seed_entropy.tsv", True)
    seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False),
                     prefix, "max_sp_seed_entropy.tsv", True)
    seed_set_entropy(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large, do_smems=False), prefix,
                                             "max_sp_"+str_msl+"_seed_entropy.tsv", True)
    
    render_times("Seed entropy - illumina", prefix, [
                    [("SMEMs l≥"+str_mss, "illumina_smem_seed_entropy.tsv", "lightblue")],
                    [("SMEMs l≥"+str_msl, "illumina_smem_"+str_msl+"_seed_entropy.tsv", "blue")],
                    [("Alg. 2a l≥"+str_msl+" (SMEMs)", "illumina_mmi_to_smem_seed_entropy.tsv", "pink")],
                    [("Alg. 1 l≥"+str_mss+" (MEMs)", "illumina_mem_seed_entropy.tsv", "orange")],
                    [("Alg. 1 l≥"+str_msl+" (MEMs)", "illumina_mem_seed_entropy_"+str_msl+".tsv", "yellow")],
                    [("10,19-minimizer l=19", "illumina_minimizer_seed_entropy.tsv", "red")],
                    [("Alg. 2b l≥"+str_msl+" (max. spanning)", "illumina_mini_to_max_sp_seed_entropy.tsv", "purple")],
                    [("maximally spanning l≥"+str_mss, "illumina_max_sp_seed_entropy.tsv", "lightgreen")],
                    [("maximally spanning l≥"+str_msl+"", "illumina_max_sp_"+str_msl+"_seed_entropy.tsv", "green")],
                ], 
                "illumina_seed_entropy",
                yAxisKey="entropy",
                #yAxisKey2="percent read covered on ref",
                y_axis_log=False )
    render_times("Seed entropy - PacBio", prefix, [
                        [("SMEMs l≥19", "pacb_smem_seed_entropy.tsv", "lightblue")],
                        [("SMEMs l≥"+str_msl, "pacb_smem_"+str_msl+"_seed_entropy.tsv", "blue")],
                        [("Alg. 2a l≥"+str_msl+" (SMEMs)", "pacb_mmi_to_smem_seed_entropy.tsv", "pink")],
                        [("Alg. 1 l≥19 (MEMs)", "pacb_mem_seed_entropy.tsv", "orange")],
                        [("Alg. 1 l≥"+str_msl+" (MEMs)", "pacb_mem_seed_entropy_"+str_msl+".tsv", "yellow")],
                        [("10,19-minimizer l=19", "pacb_minimizer_seed_entropy.tsv", "red")],
                        [("Alg. 2b l≥"+str_msl+" (max. spanning)", "pacb_mini_to_max_sp_seed_entropy.tsv", "purple")],
                        [("maximally spanning l≥"+str_mss, "pacb_max_sp_seed_entropy.tsv", "lightgreen")],
                        [("maximally spanning l≥"+str_msl, "pacb_max_sp_"+str_msl+"_seed_entropy.tsv", "green")],
                    ],
                    "pacb_seed_entropy",
                    yAxisKey="entropy",
                    #yAxisKey2="percent read covered on ref",
                    y_axis_log=False )


def runtime_analysis():
    measure_time(CreateFmdIndex(), prefix, "fm_index_construction.tsv")
    measure_time(CreateMinimizerIndex(), prefix, "minimizer_index_construction.tsv")

    measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small), prefix, "smem_computation.tsv", True)
    measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_small, do_smems=False), prefix,
                 "max_sp_seed_computation.tsv", True)
    measure_time(ComputeMaxExtendedSeeds(min_seed_length=mem_size_large), prefix,
                 "smem_"+str_msl+"_computation.tsv", True)
    measure_time(ComputeMaxExtendedSeeds(do_smems=False, min_seed_length=mem_size_large), prefix,
                 "max_sp_"+str_msl+"_seed_computation.tsv", True)
    measure_time(ComputeMinimizers(), prefix, "minimizer_computation.tsv", True)
    measure_time(LumpMinimizers(), prefix, "minimizer_lumping.tsv", True)
    measure_time(ExtendMinimizers(), prefix, "minimizer_extending.tsv", True)
    measure_time(MinimizerToSmem(), prefix, "minimizer_to_smem.tsv", True)
    measure_time(MinimizerToMaxSpan(), prefix, "minimizer_to_max_span.tsv", True)

    render_times("runtimes - illumina", prefix, [
                    [("SMEMs l≥"+str_mss, "illumina_smem_computation.tsv", "blue")],
                    [("SMEMs l≥"+str_msl, "illumina_smem_"+str_msl+"_computation.tsv", "lightblue")],
                    [("maximally spanning l≥"+str_mss, "illumina_max_sp_seed_computation.tsv", "green")],
                    [("maximally spanning l≥"+str_msl, "illumina_max_sp_"+str_msl+"_seed_computation.tsv",
                      "lightgreen")],

                    [("10,19-minimizer", "illumina_minimizer_computation.tsv", "red"),
                    ("Alg. 1 (lines 1-8)", "illumina_minimizer_lumping.tsv", "orange"),
                    ("Alg. 1 (lines 9-14)", "illumina_minimizer_extending.tsv", "yellow"),
                    ("Alg. 2a (SMEMs)", "illumina_minimizer_to_smem.tsv", "pink")],

                    [("10,19-minimizer", "illumina_minimizer_computation.tsv", "red"),
                    ("Alg. 1 (lines 1-8)", "illumina_minimizer_lumping.tsv", "orange"),
                    ("Alg. 1 (lines 9-14)", "illumina_minimizer_extending.tsv", "yellow"),
                    ("Alg. 2b (max. spanning)", "illumina_minimizer_to_max_span.tsv", "purple")],
                ],
                "illumina_times",
                divide_y_by=num_illumina_reads/1000 )
    render_times("runtimes - PacBio", prefix, [
                    [("SMEMs l≥"+str_mss, "pacb_smem_computation.tsv", "blue")],
                    [("SMEMs l≥"+str_msl, "pacb_smem_"+str_msl+"_computation.tsv", "lightblue")],
                    [("maximally spanning l≥"+str_mss, "pacb_max_sp_seed_computation.tsv", "green")],
                    [("maximally spanning l≥"+str_msl, "pacb_max_sp_"+str_msl+"_seed_computation.tsv", "lightgreen")],

                    [("10,19-minimizer", "pacb_minimizer_computation.tsv", "red"),
                    ("Alg. 1 (lines 1-8)", "pacb_minimizer_lumping.tsv", "orange"),
                    ("Alg. 1 (lines 9-14)", "pacb_minimizer_extending.tsv", "yellow"),
                    ("Alg. 2a (SMEMs)", "pacb_minimizer_to_smem.tsv", "pink")],

                    [("10,19-minimizer", "pacb_minimizer_computation.tsv", "red"),
                    ("Alg. 1 (lines 1-8)", "pacb_minimizer_lumping.tsv", "orange"),
                    ("Alg. 1 (lines 9-14)", "pacb_minimizer_extending.tsv", "yellow"),
                    ("Alg. 2b (max. spanning)", "pacb_minimizer_to_max_span.tsv", "purple")],

                ],
                "pacb_times",
                divide_y_by=num_pacb_reads/1000 )
    render_times("runtime - index creation", prefix, [
                    [("FMD-Index", "illumina_fm_index_construction.tsv", "blue")],
                    [("10,19-minimizer Index", "illumina_minimizer_index_construction.tsv", "red")],
                ],
                "index_times",
                y_axis_log=True,
                divide_y_by=0.001 )

def seed_set_diff_analysis():
    compareSeedSets(prefix, "seed_set_comparison.tsv", mem_size_small, with_single=True)
    compareSeedSets(prefix, "seed_set_comparison_"+str_msl+".tsv", mem_size_large, with_single=True)
    compareSeedSets(prefix, "seed_set_comparison_max_sp.tsv", mem_size_small, do_smems=False, with_single=True)
    compareSeedSets(prefix, "seed_set_comparison_max_sp_"+str_msl+".tsv", mem_size_large, do_smems=False,
                    with_single=True)

    render_seed_set_comp("seed set comparison - illumina", prefix, "illumina_seed_set_comparison.tsv",
                        ("unique 10,19-minimizers", "unique SMEMs l≥"+str_mss, "shared seeds l≥"+str_mss),
                         divide_y_by=num_illumina_reads)
    render_seed_set_comp("seed set comparison - pacBio", prefix, "pacb_seed_set_comparison.tsv",
                        ("unique 10,19-minimizers", "unique SMEMs l≥"+str_mss, "shared seeds l≥"+str_mss),
                         divide_y_by=num_pacb_reads)
    
    render_seed_set_comp("seed set comparison - illumina", prefix, "illumina_seed_set_comparison_"+str_msl+".tsv",
                         ("unique 10,19-minimizers", "unique SMEMs l≥"+str_msl+"", "shared seeds l≥"+str_msl), divide_y_by=num_illumina_reads)
    render_seed_set_comp("seed set comparison - pacBio", prefix, "pacb_seed_set_comparison_"+str_msl+".tsv",
                         ("unique 10,19-minimizers", "unique SMEMs l≥"+str_msl, "shared seeds l≥"+str_msl), divide_y_by=num_pacb_reads)
    
    render_seed_set_comp("seed set comparison - illumina", prefix,
                         "illumina_seed_set_comparison_max_sp.tsv",
                         ("unique 10,19-minimizers", "unique max. spanning l≥"+str_mss, "shared seeds l≥"+str_mss),
                         divide_y_by=num_illumina_reads)
    render_seed_set_comp("seed set comparison - pacBio", prefix, "pacb_seed_set_comparison_max_sp.tsv",
                         ("unique 10,19-minimizers", "unique max. spanning l≥"+str_mss, "shared seeds l≥"+str_mss),
                         divide_y_by=num_pacb_reads)
    
    render_seed_set_comp("seed set comparison - illumina", prefix,
                         "illumina_seed_set_comparison_max_sp_"+str_msl+".tsv",
                         ("unique 10,19-minimizers", "unique max. spanning l≥"+str_msl, "shared seeds l≥"+str_msl),
                         divide_y_by=num_illumina_reads)
    render_seed_set_comp("seed set comparison - pacBio", prefix,
                         "pacb_seed_set_comparison_max_sp_"+str_msl+".tsv",
                         ("unique 10,19-minimizers", "unique max. spanning l≥"+str_msl, "shared seeds l≥"+str_msl),
                         divide_y_by=num_pacb_reads)


read_generation()
runtime_analysis()
seed_entropy_analysis()
seed_set_diff_analysis()