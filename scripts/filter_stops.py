import sys
import pysam


read_ids = []
for line in sys.stdin:
    C = line.split('\t')
    mybam = pysam.AlignmentFile(sys.argv[1], "rb")

    if len(C[4]) == 1:
        for pileupcolumn in mybam.pileup('Human', int(C[1])-1,int(C[1])):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del:
                    if pileupcolumn.pos == int(C[1])-1:
                        if pileupread.alignment.query_sequence[pileupread.query_position] == C[4]:
                            read_ids.append(pileupread.alignment.query_name)
    elif len(C[4]) > 1:
        for x, char in enumerate(C[4]):
            for pileupcolumn in mybam.pileup('Human', (int(C[1])-1) + x ,(int(C[1]))+ x):
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del:
                        if pileupcolumn.pos == int(C[1])-1:
                            if pileupread.alignment.query_sequence[pileupread.query_position] == C[4][x]:
                                read_ids.append(pileupread.alignment.query_name)

c = 0
P = False
with open(sys.argv[3], 'wt') as w:
    with open(sys.argv[2], 'rt') as k:
        for l in k:
            c += 1
            l = l.strip('\n')
            if c == 1:
                idx = l.split()[0][1:]
                if idx not in read_ids:
                    P = True
                    print(l)
                else:
                    w.write(l)
                    w.write('\n')
            else:
                if P:
                    print(l)
                else:
                    w.write(l)
                    w.write('\n')
            if c >= 4:
                c = 0
                P = False
