## Arabidopsis
# get overlap peas of 2reps use idr
idr --samples SRR13189582_sub_peaks.narrowPeak SRR13189581_sub_peaks.narrowPeak \
--input-file-type narrowPeak --rank signal.value --idr-threshold 0.01 \
--output-file At_idr_g1.txt

idr --samples SRR13189569_sub_peaks.narrowPeak SRR13189564_sub_peaks.narrowPeak \
--input-file-type narrowPeak --rank signal.value --idr-threshold 0.01 \
--output-file At_idr_g2.txt

awk '{OFS="\t"}{if (NR==FNR) {id=$1$2;a[id]=$4;next} \
else {id=$1$13; pos1=int(($2+$3)/2-75); pos2=pos1+150; print $1"\t"pos1"\t"pos2"\t"a[id]"\t"($7)/2}}' \
SRR13189582_sub_peaks.narrowPeak At_idr_g1.txt > At_merged_G1_sp.bed

awk '{OFS="\t"}{if (NR==FNR) {id=$1$2;a[id]=$4;next} \
else {id=$1$13; pos1=int(($2+$3)/2-75); pos2=pos1+150; print $1"\t"pos1"\t"pos2"\t"a[id]"\t"($7)/2}}' \
SRR13189569_sub_peaks.narrowPeak At_idr_g2.txt > At_merged_G2_sp.bed

# get union peaks
bedtools intersect -a At_merged_G1_sp.bed -b At_merged_G2_sp.bed -wa -u > At_G1G2_overlap.bed
bedtools intersect -a At_merged_G1_sp.bed -b At_merged_G2_sp.bed -v > At_G1_only.bed
bedtools intersect -a At_merged_G2_sp.bed -b At_merged_G1_sp.bed -v > At_G2_only.bed
cat At_G1G2_overlap.bed At_G1_only.bed At_G2_only.bed | sort | uniq > At_G1G2_union.bed

# extend to 500bp
cat At_G1G2_union.bed | awk '{OFS="\t"}{{TSS=$4; pos1=$2-175; pos2=$3+175;}{priNb $1"\t"pos1"\t"pos2"\t"$4"\t"$5} }' > At_G1G2_all.bed_500bp.bed

# get coverage of GLK ChIP-seq
bedtools coverage -a sub_peak/At_G1G2_all.bed_500bp.bed -b bam/AtGLK1_rep1_rmdup.bam > ReadsCount_AtGLK1_rep1.At_G1G2_union_500bp.txt &
bedtools coverage -a sub_peak/At_G1G2_all.bed_500bp.bed -b bam/AtGLK1_rep2_rmdup.bam > ReadsCount_AtGLK1_rep2.At_G1G2_union_500bp.txt &
bedtools coverage -a sub_peak/At_G1G2_all.bed_500bp.bed -b bam/AtGLK2_rep1_rmdup.bam > ReadsCount_AtGLK2_rep1.At_G1G2_union_500bp.txt &
bedtools coverage -a sub_peak/At_G1G2_all.bed_500bp.bed -b bam/AtGLK2_rep2_rmdup.bam > ReadsCount_AtGLK2_rep2.At_G1G2_union_500bp.txt &

cut -f1-3,6 ReadsCount_AtGLK1_rep1.At_G1G2_union_500bp.txt > ReadsCount_AtGLK1_rep1.At_G1G2_union_500bp.RC.txt
cut -f1-3,6 ReadsCount_AtGLK1_rep2.At_G1G2_union_500bp.txt > ReadsCount_AtGLK1_rep2.At_G1G2_union_500bp.RC.txt
cut -f1-3,6 ReadsCount_AtGLK2_rep1.At_G1G2_union_500bp.txt > ReadsCount_AtGLK2_rep1.At_G1G2_union_500bp.RC.txt
cut -f1-3,6 ReadsCount_AtGLK2_rep2.At_G1G2_union_500bp.txt > ReadsCount_AtGLK2_rep2.At_G1G2_union_500bp.RC.txt

paste ReadsCount_AtGLK1_rep1.At_G1G2_union_500bp.RC.txt \
ReadsCount_AtGLK1_rep2.At_G1G2_union_500bp.RC.txt \
ReadsCount_AtGLK2_rep1.At_G1G2_union_500bp.RC.txt \
ReadsCount_AtGLK2_rep2.At_G1G2_union_500bp.RC.txt > Readscount_AtGLK_500bp.txt
