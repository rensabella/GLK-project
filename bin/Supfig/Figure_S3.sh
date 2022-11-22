## Arabidopsis
# extend to 500bp
cat At_G1G2_all_peaks.bed | awk '{OFS="\t"}{{TSS=$4; pos1=$2-175; pos2=$3+175;}{priNb $1"\t"pos1"\t"pos2"\t"$4"\t"$5} }' > At_G1G2_all.bed_500bp.bed

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
