import os
import sys
# import pysam
# from pysam.libcalignmentfile import IteratorColumnRegion
#
# def pileup_truncated(bam,contig, start, stop):
#     """
#     Obtain Pysam columns only at selected region
#     """
#     has_coord, rtid, rstart, rstop = bam.parse_region(contig, start, stop)
#     yield from IteratorColumnRegion(bam,
#                                     tid=rtid,
#                                     start=rstart,
#                                     stop=rstop,truncate=True)
#
#
#
# test_reads = ["68089efe1d51-fdde-42f4-95e9-ef44e1c25c3e", "5a79a6d8-b113-460e-80f2-807086663695"]
# read_ids = set()
# for line in sys.stdin:
#     C = line.split('\t')
#     mybam = pysam.AlignmentFile(sys.argv[1], "rb")
#
#
#     if len(C[4]) == 1:
#         for pileupcolumn in pileup_truncated(mybam, 'Human', 6020, 6040):
#             for pileupread in pileupcolumn.pileups:
#                 if pileupread.alignment.query_name in test_reads:
#                     # if not pileupread.is_del and not pileupread.is_refskip:
#                         # if pileupcolumn.pos in [int(C[1])-3, int(C[1])-2, int(C[1])-1, int(C[1]), int(C[1])+1, int(C[1])+2, int(C[1])+3]:
#                         # # sys.stderr.write("readID: {}\n".format(pileupread.alignment.query_name))
#                         # sys.stderr.write("POS: {} indel: {}, is_del: {}, is_refskip: {}, query_pos: {}, query_pos_or_next: {}\n".format(pileupcolumn.pos, pileupread.indel, pileupread.is_del,  pileupread.is_head, pileupread.is_refskip, pileupread.is_tail, pileupread.level, pileupread.query_position, pileupread.query_position_or_next))
#                         # # sys.stderr.write("Change in VCF: {} -> {}\n".format(C[3], C[4]))
#                         # try:
#                         #     sys.stderr.write("{}_{}_{}\n".format(pileupcolumn.pos, pileupread.alignment.query_sequence[pileupread.query_position]))
#                         #     sys.stderr.write("Base in read is: {}\n".format(pileupread.alignment.query_sequence[pileupread.query_position]))
#                         # except:
#                         #     sys.stderr.write("Base in read is: {}\n".format("Not in read"))
#
#
#
#                         # sys.stderr.write("{} - before\n".format(pileupread.alignment.query_name))
#                         # if not pileupread.is_del:
#                         #     # sys.stderr.write("{} - is not del\n".format(pileupread.alignment.query_name))
#                         # sys.stderr.write("{} - is in column\n".format(pileupread.alignment.query_name))
#                         # sys.stderr.write("{} first block pos hit for {}\n".format(pileupread.alignment.query_name, C[1]))
#                         sys.stderr.write("Variant {}_{}_{}\n".format(C[1], C[3], C[4]))
#                         try:
#                             sys.stderr.write("{}_{}\n".format(pileupcolumn.pos, pileupread.alignment.query_sequence[pileupread.query_position]))
#                         except:
#                             sys.stderr.write("NEXT: {}_{}\n".format(pileupcolumn.pos, pileupread.alignment.query_sequence[pileupread.query_position_or_next]))
#
#                         if pileupcolumn.pos == int(C[1]):
#                             if pileupread.alignment.query_sequence[pileupread.query_position] == C[4]:
#                                 # sys.stderr.write("{} - is at base\n".format(pileupread.alignment.query_name))
#                                 # sys.stderr.write("{} first block base hit for {}\n{}\n".format(pileupread.alignment.query_name, C[1], C))
#                                 read_ids.add(pileupread.alignment.query_name)
#                     # else:
#                     #     sys.stderr.write("{} Skipped not is_del for {}\n".format(pileupread.alignment.query_name, C[1]))
#
#     elif len(C[4]) > 1:
#         for x, char in enumerate(C[4]):
#             for pileupcolumn in mybam.pileup('Human', (int(C[1])-1) + x ,(int(C[1]))+ x):
#                 for pileupread in pileupcolumn.pileups:
#                     if pileupread.alignment.query_name in test_reads:
#                         if pileupcolumn.pos in [int(C[1])-2, int(C[1])-1, int(C[1]), int(C[1])+1]:
#                             sys.stderr.write("readID: {}\n".format(pileupread.alignment.query_name))
#                             sys.stderr.write("For pos: {} type({})\n".format(pileupcolumn.pos, type(pileupcolumn.pos)))
#                             sys.stderr.write("Change in VCF: {} -> {}\n".format(C[3], C[4]))
#                             try:
#                                 sys.stderr.write("Base in read is: {}\n".format(pileupread.alignment.query_sequence[pileupread.query_position]))
#                             except:
#                                 sys.stderr.write("Base in read is: {}\n".format("Not in read"))
#     #                 # if pileupread.alignment.query_name in test_reads:
#     #                 # sys.stderr.write("{} - before\n".format(pileupread.alignment.query_name))
#     #                 if not pileupread.is_del:
#     #                     # sys.stderr.write("{} - is not del\n".format(pileupread.alignment.query_name))
#     #                     if pileupcolumn.pos == int(C[1])-1:
#     #                         # sys.stderr.write("{} - is in column\n".format(pileupread.alignment.query_name))
#     #                         # sys.stderr.write("{} second block pos hit {}\n".format(pileupread.alignment.query_name, C[1]))
#     #                         if pileupread.alignment.query_sequence[pileupread.query_position] == C[4][x]:
#     #                             # sys.stderr.write("{} - is at base\n".format(pileupread.alignment.query_name))
#     #                             # sys.stderr.write("{} second block base hit {}\n{}\n".format(pileupread.alignment.query_name, C[1], C))
#     #                             read_ids.add(pileupread.alignment.query_name)
#     #                 # else:
#     #                 #     sys.stderr.write("{} Skipped not is_del for {}\n{}\n".format(pileupread.alignment.query_name, C[1], C))

# print("with stops", read_ids)
c = 0
P = False
with open(sys.argv[2], 'wt') as w:
    with open(sys.argv[1], 'rt') as k:
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
