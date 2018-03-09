import argparse
from tools.fileOps import *
from tools.procOps import *
from glob import glob


hub_str = '''hub NOTCH2NL
shortLabel {name}
longLabel {name}
genomesFile genomes.txt
email ian.t.fiddes@gmail.com

'''

genome_str = '''genome hg38
trackDb hg38/trackDb.txt

'''

track_template = '''track {genome}
shortLabel {genome} Assembly
longLabel {genome} Assembly
compositeTrack on
dragAndDrop subTracks
visibility dense
type bam
indelDoubleInsert on
indelQueryInsert on
showNames off
bamColorMode gray
bamGrayMode aliQual

'''

bam_template = '''    track {genome}_{bin}
    shortLabel {genome} {bin} ({num_barcodes} barcodes)
    longLabel {genome} {bin} ({num_barcodes} barcodes)
    parent {genome}
    type bam
    bigDataUrl {path}

'''


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--names', nargs='+', required=True)
    parser.add_argument('--base-dir', required=True)
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--hub-name', default='NOTCH2NL Assemblies')
    return parser.parse_args()


def align_bams(path, out_dir):
    bins = glob(os.path.join(path, 'results', 'AssembleBinsStep', '*'))
    #bins = glob(os.path.join(path, '*'))
    r = {}
    for bin in bins:
        try:
            fa = glob(os.path.join(bin, 'contig.fa'))[0]
        except IndexError:
            continue
        bin = os.path.basename(bin)
        out_path = get_tmp_file(tmp_dir=out_dir, suffix='bam')
        cmd = [['bwa', 'mem', '/hive/groups/recon/10x_genomics/references/refdata-hsapiens-hg38/fasta/genome.fa', fa],
               ['samtools', 'view', '-b', '-'],
               ['sambamba', 'sort', '-o', out_path, '/dev/stdin']]
        run_proc(cmd)
        num_barcodes = os.path.join(path, 'working', 'PhaseBarcodesStep', 'bins', bin[1:] + '.bin.txt')
        #num_barcodes = glob(os.path.join(path, bin, '*.txt'))[0]
        num_barcodes = len(open(num_barcodes).readlines())
        r[bin] = [os.path.basename(out_path), num_barcodes]
    return r


def main():
    args = parse_args()
    hg38_dir = os.path.join(args.out_dir, 'hg38')
    ensure_dir(hg38_dir)

    with open(os.path.join(args.out_dir, 'hub.txt'), 'w') as outf:
        outf.write(hub_str.format(name=args.hub_name))
    with open(os.path.join(args.out_dir, 'genomes.txt'), 'w') as outf:
        outf.write(genome_str)
    with open(os.path.join(hg38_dir, 'trackDb.txt'), 'w') as outf:
        for genome in args.names:
            outf.write(track_template.format(genome=genome))
            bams = align_bams(os.path.join(args.base_dir, genome), hg38_dir)
            for bin, (path, num_barcodes) in bams.iteritems():
                outf.write(bam_template.format(genome=genome, bin=bin, path=path, num_barcodes=num_barcodes))


if __name__ == '__main__':
    main()
