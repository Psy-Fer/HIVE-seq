import os
import sys
import argparse
import traceback
from datetime import datetime


# python3 ${WORK_DIR}/scripts/filter_stops.py -s STOP_GAINS_THAT_WE_ARE_LOOKING_FOR.txt
# -v single_read.varscan_calls.tsv -f ${IN_FASTQ} -w ${WITH_STOP_FASTQ} -n ${NO_STOP_FASTQ}

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        date_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        sys.stderr.write('[{}] - error: {}\n'.format(date_time, message))
        self.print_help()
        sys.exit(2)

def print_verbose(message, level):
    '''verbose printing'''
    if VERBOSE >= level:
        date_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        sys.stderr.write('[{}] - info: {}\n'.format(date_time, message))

def print_err(message):
    '''error printing'''
    date_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    sys.stderr.write('[{}] - error: {}\n'.format(date_time, message))

def main():
    '''
    do the thing
    '''
    NAME = "filter stop_gains"
    VERSION = "0.0.1"
    global VERBOSE
    VERBOSE=1

    parser = MyParser(
        description="filter stop_gains if the variant present")
    parser.add_argument("-s", "--stops",
                        help="stops file")
    parser.add_argument("-v", "--vars",
                        help="variants called")
    parser.add_argument("-r", "--readID",
                        help="readID")
    parser.add_argument("-f", "--filter",
                        help="write readID to filter")
    # parser.add_argument("-f", "--fastq",
    #                     help="input fatq to filter")
    # parser.add_argument("-w", "--with_stop",
    #                     help="output fastq with stops")
    # parser.add_argument("-n", "--no_stop",
    #                     help="output fastq no stops")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    stop_variant_set = set()
    with open(args.stops, 'rt') as f:
        for l in f:
            l = l.strip("\n")
            stop_variant_set.add(l)

    print_verbose(stop_variant_set, 1)

    read_variant_set = set()
    with open(args.vars, 'rt') as f:
        for l in f:
            l = l.strip("\n")
            l = l.split("\t")
            read_variant_set.add(l[0])

    INLCUDE = True
    with open(args.filter, 'a') as w:
        for v in stop_variant_set:
            if v in read_variant_set:
                INLCUDE = False
                print_verbose(v, 1)
        if INLCUDE:
            w.write(args.readID)
            w.write("\n")



if __name__ == '__main__':
    main()
