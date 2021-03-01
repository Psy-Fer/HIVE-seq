import os
import sys
import argparse
'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2019

    script description

    ----------------------------------------------------------------------------
    version 0.0 - initial



    TODO:
        -

    ----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2019 James Ferguson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
'''


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
        description="script name - script description")
    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("-s", "--sam",
                        help="read sam file")
    parser.add_argument("-l", "--long",
                        help="longest alignment file")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    readIDs = {}
    head = True
    with open(args.long, 'r') as f:
        for k in f:
            if head:
                head = False
                continue
            k = k.strip('\n')
            k = k.split('\t')
            readID = k[0]
            length = int(k[1])
            start = int(k[2])
            cigar = k[3]
            readIDs[readID] = [length, start, cigar]

    header = []
    with open(args.sam, 'r') as f:
        for k in f:
            j = k.strip('\n')
            if j[0] == "@":
                print(j)
                continue
            l = j.split('\t')
            R = l[0]
            S = int(l[3])
            C = l[5]
            diff = S - readIDs[R][1]
            if l[0] in readIDs:
                if diff < 2 and diff > -2:
                    if readIDs[R][2] in C:
                        print(j)


if __name__ == '__main__':
    main()
