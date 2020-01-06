from MA import *

## 
# @brief output folder:
# @note needs backslash at the end
# folder must exist
prefix = "/MAdata/transform_k_mers_into_smems/human_2/"

##
# @brief prefix of MA index
reference_genome_path = "/MAdata/genome/human/GRCh38.p12/ma/genome"

##
# @brief survivor and dwgsim config
dwgsim_str = "/usr/home/markus/workspace/DWGSIM/dwgsim"
num_illumina_reads = 10000
illumina_read_size = 250

survivor_str = "~/workspace/SURVIVOR/Debug/SURVIVOR simreads_n "
survivor_error_profile = "~/workspace/SURVIVOR/HG002_Pac_error_profile_bwa.txt"
num_pacb_reads = 1000

##
# @brief minimizers configuration
def get_mmi_parameter_set():
    p_m = ParameterSetManager()
    p_m.by_name("Number of Threads").set(1)
    #p_m.by_name("Minimizers - flag").set(1) # enable hcp optimization
    p_m.by_name("Minimizers - k").set(19)
    p_m.by_name("Minimizers - w").set(10)
    p_m.by_name("Use all Processor Cores").set(False)
    return p_m

##
# @brief MEM configuration
max_ambiguity_fmd = 200

mem_size_small = 19
str_mss = str(mem_size_small)

mem_size_large = 29
str_msl = str(mem_size_large)

##
# @brief x-axis configuration
start_size = 1000
stop_size = 10**7 # 3*10**9
num_steps = 10 # 100

##
# @brief save plots as svg files
# @note This requires selenium and phantomjs-prebuilt otherwise the application will crash
save_plots = False