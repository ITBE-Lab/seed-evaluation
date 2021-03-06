# A performant bridge between fixed-size and variable-size seeding
This repository contains the scripts for the experiments performed in
[A performant bridge between fixed-size and variable-size seeding](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03642-y " ").

| :exclamation: The pseudocode of Algorithm 2b contains an error: The Build-Max-Heap operation in line 7 should be "descending by *l*; for equal *l* descending by *q*" (more info below) |
|-----------------------------------------|


## Requirements

| Name | Linux install command | Recommended version | Remarks |
|------|-----------------------|---------------------|-------------|
| This repo | `git clone https://github.com/ITBE-Lab/seed-evaluation` | latest | Core evaluation scripts. |
| DWGSIM | `sudo apt-get install zlib1g-dev; sudo apt-get install libncurses5-dev; git clone --recursive https://github.com/nh13/DWGSIM; make -j$(nproc)` | 0.1.11 | Tool for generating Illumina reads. |
| SURVIVOR | `git clone https://github.com/ITBE-Lab/SURVIVOR; cd SURVIVOR/Debug; make -j$(nproc); cd ..; unzip *.zip` | 1.0.5 | Tool for generating PacBio reads. (Modified by us, to generate specific amounts of reads.) |
| Python 3 | `sudo apt-get install python3` | 3.5.3 | Python 3 environment. |
| Bokeh | `sudo apt-get install python3-pip; pip3 install bokeh` | 1.4.0 | Plotting library. |
| MA - The Modular Aligner | see below | 1.1.1-ef9ab22 | C++ library implementing all runtime critical code. |
| cmake | `sudo apt-get install cmake` | 3.13.2 | Compiling MA. |

Our testing environment: Debian GNU/Linux with a 4.9.0 kernel.

## Installing MA - The Modular Aligner

The MA github page can be found [here](https://github.com/ITBE-Lab/MA "The Modular Aligner"). \
The following sequence of commands creates the MA library:
```
git clone https://github.com/ITBE-Lab/MA
git checkout b7cf5e7            # commit used for experiments
mkdir build
cd build
cmake -DWITH_PYTHON=ON ../MA/   # with python required for evaluation scripts
make -j$(nproc)
cd ..
export PYTHONPATH=$PYTHONPATH:`pwd`/build:`pwd`/MA/python   # setup system environment
```

Type `./build/maCMD` for checking if MA was built successfully. \
If you get an error during the cmake setup or compilation, here are some things that might have gone wrong:
- MA is written in C++17, so you will need an appropriate compiler. We recommend GCC 6.3.0 or above.
- Multiple Python 3 instances on your system can confuse cmake.
- **The export command (last line in the above script) is not persistent between different terminals and logins.**

## Configuring the Python scripts

`config.py` contains the configuration:
- Set `dwgsim_str`, `survivor_str` and `survivor_error_profile` to the appropriate paths for your system.
- Set `prefix` to the folder that shall contain the temporary data and output data.
- Create an FMD-Index via `./build/maCMD --Create_Index <fasta_file_name>,<output_folder>,<index_name>`, where
`<fasta_file_name>` is the FASTA file containing the reference genome. \
Set `reference_genome_path` to `<output_folder>/<index_name>`.

## Running the experiments

Once everything is configured, you can run `python3 compute_times.py`. \
This will trigger 4 functions (very bottom of the script):
- `read_generation`: Generates the reads for all experiments and saves them in the `prefix` folder. Hence, this needs to be run only once and can be removed once the reads have been generated.
- `runtime_analysis`: Performs the time evaluation (Figure 4 of the manuscript) and generates all indices.
- `seed_entropy_analysis`: Performs the entropy analysis (Figure 5 of the manuscript).
- `seed_set_diff_analysis`: Performs the seed set difference analysis (Figure 3 of the manuscript).

Tip: Double clicking on a plot will toggle the legend visibility.

## Algorithm 2b

The published pseudocode for extracting maximal spanning seeds from MEMs contains an error: \
The Build-Max-Heap operation in line 7 should be "descending by *l*; for equal *l* descending by *q*" so that the 
first element in the max-heap is the largest seed (there can be multiple seeds with the same size) that reaches the furthest right. \
This error was only in the pseudocode; the actual implementation and measurements are correct.

Further, a max-heap is not actually required. Instead a single iteration over *T* is enough to extract all relevant seeds:
![Alt text](mems_to_max_sp.png?raw=true "Algorithm 2b")

We would like to thank Roman Cheplyaka for pointing out our mistake as well as this optimization.