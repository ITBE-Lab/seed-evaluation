from config import *
import os


def backwardMEM(reference_filename, queries_filename, min_seed_len):
    os.system(backward_mem_str + "-l " + str(min_seed_len) + " -b " + reference_filename + " " +
              queries_filename + " > /dev/null 2>&1 ")