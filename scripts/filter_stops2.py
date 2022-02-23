import os
import sys
import argparse
import traceback
from datetime import datetime


import pysam


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
                        help="input fastq to filter")
    parser.add_argument("-t", "--readIDs",
                        help="readID list")
    # parser.add_argument("-r", "--read_bases",
    #                     help="readID truncated bases")
    parser.add_argument("-b", "--bam",
                        help="bam file")
    parser.add_argument("-w", "--with_stop",
                        help="output fastq with stops")
    parser.add_argument("-n", "--no_stop",
                        help="output fastq no stops")

    args = parser.parse_args()




    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    check_positions = []
    with open(args.readIDs, 'r') as f:
        for l in f:
            l = l.strip("\n")
            # readID = l[0]
            # read_pos = l[4]
            # ref_pos = l[7]
            ref_pos = int(l)
            # if readID not in filter_set:
            #     filter_set[l[0]] = {}

            # filter_set[l[0]][ref_pos] = read_pos
            check_positions.append(ref_pos)

    print(check_positions)


    read_dic = {}
    offset = 0
    F = pysam.Samfile(args.bam, "rb")
    for r in F.fetch():
        name = r.query_name
        # if name not in filter_set:
        #     continue
        if name not in read_dic:
            read_dic[name] = {}
        seq = r.query_alignment_sequence
        start = r.query_alignment_start
        end = r.query_alignment_end
        positions = r.get_aligned_pairs() # (Q, R) tupple, could be None
        # print(positions)
        new_seq_list = []
        # print(name)
        # print("Q:", len(seq), "R:", len(positions), "q_start:", start, end)
        # print(seq)
        base_skip = False
        skip_count = 0
        for k in range(len(positions)):
            n, m = positions[k]
            if k < start:
                continue
            if k > end:
                break
            if n is None:
                base_skip = True
                skip_count += 1
                new_seq_list.append("-")

            else:
                if positions[k][1] in check_positions:
                    # print(k, skip_count)
                    # print(seq[k-skip_count-start], positions[k], "<========================", check_positions)
                    new_seq_list.append(seq[k-skip_count-start])
                else:
                    # print(k, skip_count)
                    # print(seq[k-skip_count-start], positions[k])
                    new_seq_list.append(seq[k-skip_count-start])
        new_seq = "".join(new_seq_list)

        base_skip = False
        for i in range(len(positions)):
            Q, R = positions[i] # could be None..need to account for this
            if i < start:
                continue
            if i > end:
                break
            if Q is None:
                continue
            if R in check_positions:
                motif = new_seq[(i-start)-3:(i-start)+2]
                neg_offset = 0
                pos_offset = 0
                while "-" in motif:
                    motif = list(motif)
                    # print("motif len", len(motif))
                    # print(motif)
                    for j in range(len(motif)):
                        # print("index:", j)
                        # print(motif)
                        if motif[j] == "-":
                            if j < 2:
                                neg_offset += 1
                                motif.pop(j)
                                left_add = new_seq[(i-start)-3-neg_offset]
                                # print("motif", type(motif), "left_add", type(left_add))
                                motif = list(left_add) + list(motif)
                            if j >= 2:
                                pos_offset += 1
                                motif.pop(j)
                                right_add = [new_seq[(i-start)+2+pos_offset]]
                                motif = list(motif) + list(right_add)
                    motif = "".join(motif)

                read_dic[name][R] = motif

                # print(name, "Q:", Q, "R:", R, R, start, motif)
            offset = 0

    # print(read_dic)
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
            line1 = l
        elif c == 2:
            check = False
            for p in check_positions:
                pos = p
                # if idx not in read_dic: # TESTING ONLY!!!!!!!!
                #     continue
                if pos in read_dic[idx]:
                    print(idx, pos)
                    if check_stop(read_dic[idx][pos]):
                        check = True
                        print("True")
                    else:
                        print("False")
                else:
                    print("Position not found!!!!")
            if check:
                W.write(line1)
                W.write('\n')
                W.write(l)
                W.write('\n')
            else:
                check = False
                N.write(line1)
                N.write('\n')
                N.write(l)
                N.write('\n')
        else:
            if check:
                W.write(l)
                W.write('\n')
            else:
                N.write(l)
                N.write('\n')
        if c >= 4:
            c = 0
            check = False

    F.close()
    W.close()
    N.close()

def check_stop(seq):
    '''
    check that stop is still a stop in case another variant modifies the codon
    N = position in the read sequence.
    Look at the 2 bases either side to get context for all 3bp codons
    go through steps 0, 1, 2 to get the 3 possible codons for the change.
    Check if any of them meet the stop definitions
    '''
    # motif = seq[N-2:N+3]
    # print(seq[N])
    print(seq)
    for i in [0, 1, 2]:
        # if seq[i:i+3] in ["TGA", "TAG", "TAA", "TCA", "CTA", "TTA"]:
        if seq[i:i+3] in ["TGA", "TAG", "TAA"]:
            return True
    return False

if __name__ == '__main__':
    main()
