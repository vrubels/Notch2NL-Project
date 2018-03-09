import pysam
import os
import pandas as pd
from pyfaidx import Fasta
from collections import *
from tools.intervals import *
from tools.misc import *
from tools.procOps import *
from tools.fileOps import *
from tools.bio import *
from tools.psl import *
from itertools import *
import bisect
import argparse


class MsaRecord(object):
    """Convenience holder for alignment information"""
    def __init__(self, msa_start, msa_stop, qname, seq):
        self.msa_start = msa_start
        self.msa_stop = msa_stop
        self.qname = qname
        self.seq = seq


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n2nl-aln', default='/hive/users/ifiddes/notch_brian_viz/try_bwa/notch2nl_alignment.fa')
    parser.add_argument('--hg38-bam', default='/hive/users/ifiddes/notch_brian_viz/try_bwa/hg38_mapped.bam')
    parser.add_argument('--n2-regions', default='/hive/users/ifiddes/notch_brian_viz/try_bwa/n2_regions.bed')
    parser.add_argument('--hg38-bams', required=True, nargs='+')
    parser.add_argument('--output-fasta', required=True, type=argparse.FileType('w'))
    return parser.parse_args()


def find_closest(numeric_list, query_number):
    """
    Given a list of numbers, and a single query number, find the number in the sorted list that is numerically
    closest to the query number. Uses list bisection to do so, and so should be O(log n)
    """
    sorted_numeric_list = sorted(numeric_list)
    pos = bisect.bisect_left(sorted_numeric_list, query_number)
    if pos == 0:
        return sorted_numeric_list[0]
    if pos == len(sorted_numeric_list):
        return sorted_numeric_list[-1]
    before = sorted_numeric_list[pos - 1]
    after = sorted_numeric_list[pos]
    if after - query_number < query_number - before:
        return after
    else:
        return before


def find_overlap(b1, b2, spacer=50):
    delta = b1.msa_stop - b2.msa_start
    b1_seq = b1.seq[-(delta + spacer):]
    b2_seq = b2.seq[:delta + spacer]
    aln = perform_aln(b1_seq, b2_seq)
    if len(aln) == 0:
        s = b1.seq[:-delta]
        #assert len(s) > 0, (b1.qname, b2.qname)
        return s
    else:
        psl = sorted([PslRow(x.split('\t')) for x in aln], key=lambda x: x.coverage)[-1]
        cutoff = psl.query_coordinate_to_target(psl.q_end - 1) + 1
        return b1.seq[:-(cutoff - spacer)]


def perform_aln(b1_seq, b2_seq):
    with TemporaryFilePath() as b1_f, TemporaryFilePath() as b2_f:
        with open(b1_f, 'w') as b1_f_h:
            write_fasta(b1_f_h, 'b1', b1_seq)
        with open(b2_f, 'w') as b2_f_h:
            write_fasta(b2_f_h, 'b2', b2_seq)
        cmd = ['blat', b1_f, b2_f, '-noHead', '/dev/stdout']
        return call_proc_lines(cmd)[:-2]


def construct_hg38_map(n2nl_aln, hg38_bam):
    """Constructs a map of hg38 position -> sequence alignment position -> MSA position"""
    # construct sequence alignment position -> MSA position map using the MSA
    aln_f = Fasta(n2nl_aln)
    seq_aln_map = defaultdict(dict)
    for name, seq in aln_f.iteritems():
        seq_pos = 0
        for aln_pos, x in enumerate(str(seq)):
            seq_aln_map[name][seq_pos] = aln_pos
            if x != '-':
                seq_pos += 1

    # find maximum position for reversing negative strand
    max_pos = {x: max(y.keys()) for x, y in seq_aln_map.iteritems()}

    # construct a hg38 -> sequence positions using the sequences trivially mapped back to hg38
    hg38_map = {}
    for rec in pysam.Samfile(hg38_bam):
        m = {y: x for x, y in rec.aligned_pairs}
        # invert positions for negative strand genes
        if rec.qname in ['NOTCH2', 'NOTCH2NL-A', 'NOTCH2NL-B']:
            m = {x: max_pos[rec.qname] - y for x, y in m.iteritems()}
        hg38_map[rec.qname] = m

    # construct a table mapping each alignment position to all hg38 positions
    r = defaultdict(dict)
    for name, pos_map in hg38_map.iteritems():
        for hg38_pos, seq_pos in pos_map.iteritems():
            aln_pos = seq_aln_map[name][seq_pos]
            r[name][aln_pos] = hg38_pos

    # now invert this map, so that we have our hg38 -> aln map
    final_map = {}
    for name in r:
        for aln_pos in r[name]:
            hg38_pos = r[name][aln_pos]
            assert hg38_pos not in final_map
            final_map[hg38_pos] = aln_pos

    return final_map


def load_alignments(bam, regions_of_interest, position_map, sorted_positions):
    """
    1) Load all alignments of contigs to hg38
    2) filter for those that overlap the MSA
    3) sort them by MSA positions based on the position map
    4) Filter for entirely overlapping
    5) return sorted lists of MsaRecords
    """
    msa_blocks = []
    for aln in pysam.Samfile(bam):
        if aln.is_unmapped:
            continue
        # make sure this alignment overlaps notch
        start = aln.reference_start
        end = aln.reference_end
        c = ChromosomeInterval('chr1', start, end, '.')
        # filter for overlapping our regions and not being too short
        if not interval_not_intersect_intervals(regions_of_interest, c) and len(c) > 100:
            # find the closest positions in hg38 that are present in our MSA
            closest_start = find_closest(sorted_positions, start)
            closest_stop = find_closest(sorted_positions, end)
            msa_start = position_map[closest_start]
            msa_stop = position_map[closest_stop]
            if msa_start > msa_stop:  # handle negative strand
                msa_start, msa_stop = msa_stop, msa_start
                seq = reverse_complement(aln.seq)
            else:
                seq = aln.seq
            b = MsaRecord(msa_start, msa_stop, aln.qname, seq)
            msa_blocks.append(b)

    # sort by MSA position
    sorted_msa_blocks = sorted(msa_blocks, key=lambda x: x.msa_start)

    # filter blocks for those that are entirely a subset of another
    # this removes misassemblies
    intervals = [ChromosomeInterval('', x.msa_start, x.msa_stop, '.', x) for x in sorted_msa_blocks]
    bad_intervals = set()
    for i1, i2 in combinations(intervals, 2):
        if i1.subset(i2):
            bad_intervals.add(i1)
        elif i2.subset(i1):
            bad_intervals.add(i2)
    filtered_msa_blocks = [x.data for x in intervals if x not in bad_intervals]

    return filtered_msa_blocks


def scaffold_alignments(filtered_msa_blocks):
    """
    Takes a sorted list of MsaRecord objects and scaffolds them
    """
    seq = []
    # construct a pairwise iterator over the filtered_msa_blocks
    a, b = tee(filtered_msa_blocks)
    _ = next(b, None)
    for i, (b1, b2) in enumerate(izip(a, b)):
        if b2.msa_start < b1.msa_stop:  # we have an overlap, resolve via pairwise alignment
            seq.append(find_overlap(b1, b2))
        elif b1.msa_stop == b2.msa_start:  # exactly contiguous, just add b1 and continue
            seq.append(b1.seq)
        else:  # we have a unknown gap, so add b1 then a gap
            seq.append(b1.seq)
            seq.append(''.join(['N'] * 100))
        # add the final b2 sequence
    seq.append(b2.seq)
    return ''.join(seq)


if __name__ == '__main__':
    args = parse_args()
    # load the intervals we are interested in
    start_stop_positions = {x.split()[3]: (x.split()[1], x.split()[2]) for x in open(args.n2_regions)}
    start_stop_positions = {x: map(int, y) for x, y in start_stop_positions.iteritems()}

    # load the map of hg38 positions to alignment positions
    position_map = construct_hg38_map(args.n2nl_aln, args.hg38_bam)
    # construct a sorted list of hg38 positions
    sorted_positions = sorted(position_map.keys())

    # construct our regions of interest
    regions_of_interest = []
    for start, stop in start_stop_positions.itervalues():
        regions_of_interest.append(ChromosomeInterval('chr1', start, stop, '.'))

    for bam in args.hg38_bams:
        name = os.path.basename(bam).split('.')[0]
        filtered_msa_blocks = load_alignments(bam, regions_of_interest, position_map, sorted_positions)
        seq = scaffold_alignments(filtered_msa_blocks)
        write_fasta(args.output_fasta, name, seq)
