set -e
g=$1
ploidy=$2
DIR=$3

realign='/hive/users/ifiddes/notch2nl_berkeley_data/extract_realign_fix_tags_low_memory.py --regions chr1:119989248-120190000 chr1:149328818-149471561 chr1:148600079-148801427 chr1:146149145-146328264 chr1:120705669-120801220'

targets=/hive/users/ifiddes/notch2nl_berkeley_data/new_consensus/unmasked_merged.bed
unmasked=/hive/users/ifiddes/notch2nl_berkeley_data/new_consensus/extended_alignment_consensus.fa
cnv_map=/hive/users/ifiddes/notch2nl_berkeley_data/new_consensus/normal_diploid_copies.bed

bam=${g}/outs/phased_possorted_bam.bam

if [ ! -e ${bam} ]; then
  bam=${g}/outs/possorted_bam.bam
fi

python ${realign} --consensus_ref ${unmasked} ${bam} ${g}/${g}.consensus_mapped.sorted.bam ${g}/${g}.fastq --paired

freebayes -f ${unmasked} --ploidy ${ploidy} --min-alternate-fraction 0.05 -k -j --min-coverage 50 -i -u -0 --cnv-map ${cnv_map} \
${g}/${g}.consensus_mapped.sorted.bam > ${g}/${g}.${ploidy}.variant_calls.vcf

python /hive/users/ifiddes/notch2nl_berkeley_data/sortfq.py ${g}/${g}.fastq ${g}/${g}.sorted.fastq

#cd ${g}
athena ${g}/full/config.json
