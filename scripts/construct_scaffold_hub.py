import argparse
from tools.fileOps import *
from tools.procOps import *
from tools.bio import *
from pyfasta import Fasta


hub_str = '''hub {0}
shortLabel {1}
longLabel {1}
genomesFile genomes.txt
email NoEmail

'''


genomes_str = '''genome Notch2NL_consensus
twoBitPath consensus/consensus.2bit
trackDb consensus/trackDb.txt
organism Notch2NL_consensus
description Notch2NL_consensus
scientificName Notch2NL_consensus
defaultPos Notch2NL_consensus:1-1000000

'''


track_str = '''track {0}
shortLabel {0} alignment
longLabel {0} alignment
type bigPsl
baseColorUseSequence lfExtra
baseColorDefault diffBases
showDiffBasesAllScales on
bigDataUrl {0}.bb
visibility pack
priority {1}

'''

annot_str = '''track notch2nl
shortLabel NOTCH2NL
longLabel NOTCH2NL
bigDataUrl annotations.bb
type bigBed 12
visibility pack
priority 1

'''


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scaffolded', required=True)
    parser.add_argument('--out-hub', required=True)
    parser.add_argument('--consensus', default='/hive/users/ifiddes/notch_brian_viz/try_lastz/consensus.fa')
    parser.add_argument('--consensus_2bit', default='/hive/users/ifiddes/notch_brian_viz/try_lastz/consensus.2bit')
    parser.add_argument('--consensus_sizes', default='/hive/users/ifiddes/notch_brian_viz/try_lastz/consensus.chrom.sizes')
    parser.add_argument('--big-psl-as', default='/cluster/home/ifiddes/kent/src/hg/lib/bigPsl.as')
    parser.add_argument('--hub-name', default='NOTCH2NL consensus mapping')
    parser.add_argument('--annotation-bb', default='/hive/users/ifiddes/notch_brian_viz/try_lastz/hub/consensus/annotations.bb')
    return parser.parse_args()


def construct_big_psl(name, seq, consensus_fa, consensus_2bit, consensus_sizes, big_psl_as, out_big_psl):
    """Constructs a bigPsl mapping this sequence to the consensus using lastz and chains"""
    with TemporaryFilePath() as seq_sizes, TemporaryFilePath() as seq_2bit, TemporaryFilePath() as tmp_psl, \
        TemporaryFilePath() as seq_fa:
        write_fasta(seq_fa, name, str(seq))
        cmd = ['faToTwoBit', seq_fa, seq_2bit]
        run_proc(cmd)
        cmd = ['twoBitInfo', seq_2bit, seq_sizes]
        run_proc(cmd)
        cmd = [['lastz', consensus_fa, seq_fa, '--strand=plus', '--output=/dev/stdout'],
               ['lavToPsl', '/dev/stdin', '/dev/stdout'],
               ['axtChain', '-linearGap=medium', '-psl', '/dev/stdin', consensus_2bit, seq_2bit, '/dev/stdout'],
               ['chainPreNet', '/dev/stdin', consensus_sizes, seq_sizes, '/dev/stdout'],
               ['chainToPsl', '/dev/stdin', consensus_sizes, seq_sizes, consensus_2bit, seq_2bit, '/dev/stdout'],
               ['pslFilter', '/dev/stdin', tmp_psl]]
        run_proc(cmd)
        # now we need 1 last temp file
        with TemporaryFilePath() as tmp_bgp:
            cmd = [['pslToBigPsl', tmp_psl, 'stdout', '-fa={}'.format(seq_fa)],
                   ['sort', '-k1,1', '-k2,2n']]
            run_proc(cmd, stdout=tmp_bgp)
            cmd = ['bedToBigBed', '-type=bed12+13', '-tab', '-as={}'.format(big_psl_as), tmp_bgp, consensus_sizes,
                   out_big_psl]
            run_proc(cmd)


if __name__ == '__main__':
    args = parse_args()
    ensure_dir(os.path.join(args.out_hub, 'consensus'))

    with open(os.path.join(args.out_hub, 'genomes.txt'), 'w') as outf:
        outf.write(genomes_str)

    with open(os.path.join(args.out_hub, 'hub.txt'), 'w') as outf:
        outf.write(hub_str.format('_'.join(args.hub_name.split(' ')), args.hub_name))

    with open(os.path.join(args.out_hub, 'consensus', 'trackDb.txt'), 'w') as outf:
        outf.write(annot_str)
        os.link(args.annotation_bb, os.path.join(args.out_hub, 'consensus/annotations.bb'))
        os.link(args.consensus_2bit, os.path.join(args.out_hub, 'consensus/consensus.2bit'))
        f = Fasta(args.scaffolded)
        for i, (name, seq) in enumerate(sorted(f.iteritems()), 2):
            out_big_psl = os.path.join(args.out_hub, 'consensus', '{}.bb'.format(name))
            big_psl = construct_big_psl(name, seq, args.consensus, args.consensus_2bit, args.consensus_sizes,
                                        args.big_psl_as, out_big_psl)
            outf.write(track_str.format(name, i))
