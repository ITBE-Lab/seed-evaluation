from MA import *

## 
# @brief output folder:
# @note needs backslash at the end
# folder must exist
prefix = "/MAdata/transform_k_mers_into_smems/human_ccs/"

##
# @brief prefix of MA index
reference_genome_path = "/MAdata/genome/human/GRCh38.p12/ma/genome"
reference_genome_fasta = "/MAdata/genome/human/GRCh38.p12/fasta/genome.fna"

##
# @brief survivor and dwgsim config
dwgsim_str = "/usr/home/markus/workspace/DWGSIM/dwgsim"
num_illumina_reads = 10000
illumina_read_size = 250

survivor_str = "~/workspace/SURVIVOR/Debug/SURVIVOR simreads_n "
survivor_error_profile = "~/workspace/SURVIVOR/HG002_Pac_error_profile_bwa.txt"
num_pacb_reads = 1000

##
# 
backward_mem_str = "~/workspace/backwardMEM/backwardMEM/backwardMEM8 "

##
# @brief minimizers configuration
def get_mmi_parameter_set():
    p_m = ParameterSetManager()
    p_m.by_name("Number of Threads").set(1)
    #p_m.by_name("Minimizers - flag").set(1) # enable hcp optimization
    p_m.by_name("Minimizer Size").set(19) # k
    p_m.by_name("Minimizer Window Size").set(10) # w
    p_m.by_name("Use all Processor Cores").set(False)
    return p_m

##
# @brief occurence filter configuration configuration
max_ambiguity_fmd = 2000

mem_size_small = 19 # == k
str_mss = str(mem_size_small)

mem_size_large = 28 # == w + k - 1
str_msl = str(mem_size_large)


x_axis_unit = "read_noise" # or "genome_section_size" (to run the old analysis)


if x_axis_unit == "genome_section_size":
    ##
    # @brief x-axis configuration old analysis
    start_size = 1000
    stop_size = 3*10**9
    num_steps = 100
if x_axis_unit == "read_noise":
    ##
    # @brief x-axis configuration new analysis
    start_size = 0
    stop_size = 1.1
    num_steps = 10

##
# @brief save plots as svg files
# @note This requires selenium and phantomjs-prebuilt otherwise the application will crash
save_plots = False