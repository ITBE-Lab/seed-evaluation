# A performant bridge between fixed-size and variable-size seeding
This repository contains the experiments performed in
[A performant bridge between fixed-size and variable-size seeding](https://doi.org/10.1101/825927 "bioRxiv - preprint").

## Requirements

| Name | Linux install command | Recommended version | Remarks |
|------|-----------------------|---------------------|-------------|
| This repo | `git clone https://github.com/ITBE-Lab/seed-evaluation` | latest | Core evaluation scripts. |
| DWGSIM | `sudo apt-get install zlib1g-dev; git clone --recursive https://github.com/nh13/DWGSIM; make -j32` | 0.1.11 | Tool for generating Illumina reads. |
| SURVIVOR | `git clone https://github.com/ITBE-Lab/SURVIVOR; cd SURVIVOR/Debug; make -j32; cd ..; unzip *.zip` | 1.0.5 | Tool for generating PacBio reads. (Modified by us, to generate specific amounts of reads.) |
| Python 3 | `sudo apt-get install python3` | 3.5.3 | Python 3 environment. |
| Bokeh | `sudo apt-get install python3-pip; pip3 install bokeh` | 1.4.0 | Plotting library. |
| MA - The Modular Aligner | see below | 1.1.1-1f5a63f | C++ library needed by core evaluation scripts. |

Our testing environment: Ubuntu 18.04.3 LTS.

## Installing MA - The Modular Aligner

The MA github page can be found [here](https://github.com/ITBE-Lab/MA "The Modular Aligner"). \
The following sequence of commands creates the MA library:
```
git clone https://github.com/ITBE-Lab/MA
git checkout 1f5a63f            # commit used for experiments
mkdir build
cd build
cmake -DWITH_PYTHON=ON ../MA/   # with python required for evaluation scripts
make -j32
cd ..
export PYTHONPATH=$PYTHONPATH:`pwd`/build:`pwd`/MA/python   # setup system environment
```

Type `./build/maCMD` for checking if MA was built successfully. \
If you get an error during the cmake setup or compilation, here are some things that might have gone wrong:
- MA is written in C++17, so you will need an appropriate compiler. We recommend GCC 6.3.0 or above.
- Multiple Python 3 instances on your system can confuse cmake.
- The export command (last line in the above script) is not persistent between different terminals and logins.

## Configuring the python scripts

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

Tip: Double clicking on a plot will toggle its legend.