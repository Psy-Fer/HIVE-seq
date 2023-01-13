import os
import sys
import argparse
import traceback
from datetime import datetime

# python3 ${WORK_DIR}/scripts/filter_stops.py -f ${IN_FASTQ} -t ${FILER_READS} -w ${WITH_STOP_FASTQ} -n ${NO_STOP_FASTQ}

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

    parser = MyParser(
        description="filter stop_gains if the variant present")
    parser.add_argument("-f", "--fastq",
                        help="input fatq to filter")
    parser.add_argument("-t", "--readIDs",
                        help="readID list")
    parser.add_argument("-w", "--with_stop",
                        help="output fastq with stops")
    parser.add_argument("-n", "--no_stop",
                        help="output fastq no stops")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    filter_set = set()
    with open(args.readIDs, 'r') as f:
        for l in f:
            l = l.strip("\n")
            l = l.split("\t")
            filter_set.add(l[0])


    F = open(args.fastq, 'r')
    W = open(args.with_stop, 'w')
    N = open(args.no_stop, 'w')
    c = 0
    P = False
    for l in F:
        c += 1
        l = l.strip('\n')
        if c == 1:
            idx = l.split()[0][1:]
            if idx in filter_set:
                P = True
                W.write(l)
                W.write('\n')
            else:
                N.write(l)
                N.write('\n')
        else:
            if P:
                W.write(l)
                W.write('\n')
            else:
                N.write(l)
                N.write('\n')
        if c >= 4:
            c = 0
            P = False

    F.close()
    W.close()
    N.close()

if __name__ == '__main__':
    main()
