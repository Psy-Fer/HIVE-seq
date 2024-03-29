import os
import sys
import argparse
import numpy as np
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
    parser.add_argument("-i", "--input",
                        help="file to read")
    parser.add_argument("-o", "--output",
                        help="output file")


    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    '''
    For each row if column 4 value is less than 10 = TRUE enter the value in column 3
    If column 34>column 40 > column 46 > column 52 = TRUE enter A
    If column 40 > column 34 & 46 & 52 = TRUE enter C
    If column 46 > column 34 & 40 & 52 = TRUE enter T
    If column 52 > column 34 & 40 & 46 = TRUE enter G
    '''
    opt = ['A', 'C', 'T', 'G']
    cons = []
    head = True
    with open(args.input, 'r') as f:
        for l in f:
            if head:
                head = False
                continue
            l = l.strip('\n')
            l = l.split('\t')
            if int(l[3]) < 10:
                cons.append(l[2])
            else:
                # A, C, T, G
                scores = [int(l[33]), int(l[39]), int(l[45]), int(l[51])]
                max = np.argmax(scores)
                base = opt[max]
                cons.append(base)

    C = 1
    with open(args.output, 'w') as w:
        w.write('>Human\n')
        for i in cons:
            w.write(i)
            if C >= 60:
                w.write('\n')
                C = 0
            C += 1


if __name__ == '__main__':
    main()
