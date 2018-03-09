#list of normals we are going to evaluate
bams = {'NA19240': '/hive/users/ifiddes/notch2nl_berkeley_data/NA19240/outs/phased_possorted_bam.bam',
       'H9': '/hive/users/ifiddes/notch2nl_berkeley_data/E2del19N_E2del68_combined_longranger/E2del68_E2del19N_combined/outs/phased_possorted_bam.bam',
       'NA12877': '/hive/users/ifiddes/notch2nl_berkeley_data/NA12877/outs/phased_possorted_bam.bam'}

#!mkdir cnv_map
import os
os.chdir('cnv_map')

# for each normal, construct a bigWig over the unique region
from tools.procOps import *
from tools.fileOps import *
unique_region = 'chr1:142785299-150598866'
size = 150598866 - 142785299
for g, bam in bams.iteritems():
    wig = g + '.bw'
    cmd = ['bamCoverage', '-b', bam, '-o', wig, '--numberOfProcessors', '20',
          '--normalizeTo1x', size, '--region', unique_region]
    run_proc(cmd)

# bring in the seg dup depth bed, cutting it up
#!grep ^chr /hive/users/ifiddes/notch2nl_berkeley_data/seg_dup_depth/chr1_seg_dup_depth_with_color.bed | cut -d $'\t' -f 1,2,3 > seg_dupes.bed
# also include the depth version for later
#!grep ^chr /hive/users/ifiddes/notch2nl_berkeley_data/seg_dup_depth/chr1_seg_dup_depth_with_color.bed | tawk '{print $1, $2, $3, $4 + 1}' > seg_dupes_with_depth.bed
# construct a bed of the complement of the seg dup depth map to get the single copy regions
#!bedtools complement -g /cluster/data/hg38/chrom.sizes -i seg_dupes.bed > seg_dup_complement.bed
# add the value to this bed
#!cat seg_dup_complement.bed | tawk '{print $1, $2, $3, 1}' > seg_dup_complement_1.bed
# combine and sort
#!cat seg_dup_complement_1.bed seg_dupes_with_depth.bed > combined.bed
#!bedSort combined.bed combined.sorted.bed

# construct the average
cmd = ['wiggletools', 'mean'] + ['{}.bw'.format(x) for x in bams]
run_proc(cmd, stdout='averaged.wig')

# scale based on seg dup depth track
cmd=['wiggletools', 'ratio', 'averaged.wig', 'combined.sorted.bed']
run_proc(cmd, stdout='averaged_scaled.wig')

# convert to bigWig for later browser usage
cmd=['wigToBigWig', 'averaged_scaled.wig', '/cluster/data/hg38/chrom.sizes', 'averaged_scaled.bw']
run_proc(cmd)

# test against inputs
for g in bams:
    cmd = ['wiggletools', 'ratio', '{}.bw'.format(g), 'averaged_scaled.bw']
    run_proc(cmd, stdout='{}_ratio.wig'.format(g))
    cmd = ['wigToBigWig', '{}_ratio.wig'.format(g), '/cluster/data/hg38/chrom.sizes', '{}_ratio.bw'.format(g)]
    run_proc(cmd)

# this looked pretty good, now lets try the patients
#patients = {'PatientA': '/hive/users/ifiddes/notch2nl_berkeley_data/PatientA/outs/phased_possorted_bam.bam',
#           'PatientB': '/hive/users/ifiddes/notch2nl_berkeley_data/PatientB/outs/phased_possorted_bam.bam',
#           'E5del42': '/hive/users/ifiddes/notch2nl_berkeley_data/E5del42/outs/phased_possorted_bam.bam',
#           'E5del61': '/hive/users/ifiddes/notch2nl_berkeley_data/E5del61/outs/phased_possorted_bam.bam',
#           'E5del54': '/hive/users/ifiddes/notch2nl_berkeley_data/E5del54/outs/phased_possorted_bam.bam',
#           'E5del58': '/hive/users/ifiddes/notch2nl_berkeley_data/E5del58/outs/phased_possorted_bam.bam'}
patients = { 'CHM1' : '/hive/users/cbosworth/170714_10X/CHM1/outs/phased_possorted_bam.bam',
            'CHM13' : '/hive/users/cbosworth/170714_10X/CHM13/outs/phased_possorted_bam.bam',
         'E2del19N' : '/hive/users/cbosworth/170714_10X/E2del19N/outs/phased_possorted_bam.bam',
          'E2del68' : '/hive/users/cbosworth/170714_10X/E2del68/outs/phased_possorted_bam.bam',
          'E2del70' : '/hive/users/cbosworth/170714_10X/E2del70/outs/phased_possorted_bam.bam',
          'E5del42' : '/hive/users/cbosworth/170714_10X/E5del42/outs/phased_possorted_bam.bam',
            'ESC12' : '/hive/users/cbosworth/170714_10X/ESC12/outs/phased_possorted_bam.bam',
          'NA12877' : '/hive/users/cbosworth/170714_10X/NA12877/outs/phased_possorted_bam.bam',
          'NA19238' : '/hive/users/cbosworth/170714_10X/NA19238/outs/phased_possorted_bam.bam',
          'NA19239' : '/hive/users/cbosworth/170714_10X/NA19239/outs/phased_possorted_bam.bam',
          'NA19240' : '/hive/users/cbosworth/170714_10X/NA19240/outs/phased_possorted_bam.bam',
         'PatientA' : '/hive/users/cbosworth/170714_10X/PatientA/outs/phased_possorted_bam.bam',
         'PatientB' : '/hive/users/cbosworth/170714_10X/PatientB/outs/phased_possorted_bam.bam',
            'SV721' : '/hive/users/cbosworth/170714_10X/SV721/outs/phased_possorted_bam.bam',
            'SV735' : '/hive/users/cbosworth/170714_10X/SV735/outs/phased_possorted_bam.bam',
            'SV780' : '/hive/users/cbosworth/170714_10X/SV780/outs/phased_possorted_bam.bam',
            'SV788' : '/hive/users/cbosworth/170714_10X/SV788/outs/phased_possorted_bam.bam',
            'SV877' : '/hive/users/cbosworth/170714_10X/SV877/outs/phased_possorted_bam.bam'}
for g, bam in patients.iteritems():
    wig = g + '.bw'
    cmd = ['bamCoverage', '-b', bam, '-o', wig, '--numberOfProcessors', '20',
          '--normalizeTo1x', size, '--region', unique_region]
    run_proc(cmd)
    cmd = ['wiggletools', 'ratio', '{}.bw'.format(g), 'averaged_scaled.bw']
    run_proc(cmd, stdout='{}_ratio.wig'.format(g))
    cmd = ['wigToBigWig', '{}_ratio.wig'.format(g), '/cluster/data/hg38/chrom.sizes', '{}_ratio.bw'.format(g)]
    run_proc(cmd)


