import os
import sys
import argparse


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def main():
    '''
    do the thing
    '''
    parser = MyParser(
        description="qfilter - filter fastq for q score usig seq_sum")
    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("-f", "--fastq",
                        help="fastq file to filter")
    parser.add_argument("-s", "--seq_sum",
                        help="seq sum file input")
    parser.add_argument("-q", "--qual", type=float, default=7.0,
                        help="quality score >= to [Default=7.0]")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    read_set = set()
    head = True
    with open(args.seq_sum, 'r') as f:
        for l in f:
            l = l.split('\t')
            if head:
                head = False
                col = False
                for i in range(len(l)):
                    if l[i] == "mean_qscore_template":
                        col = i
                if not col:
                    sys.stderr.write("column mean_qscore_template not found")
                    sys.exit()
                continue
            if float(l[col]) >= args.qual:
                read_set.add(l[1])


    c = 0
    P = False
    with open(args.fastq, 'r') as f:
        for ln in f:
            c += 1
            l = ln.strip('\n')
            if c == 1:
                idx = l.split()[0][1:]
                if idx in read_set:
                    P = True
                    print(l)
            else:
                if P:
                    print(l)
            if c >= 4:
                c = 0
                P = False



if __name__ == '__main__':
    main()
