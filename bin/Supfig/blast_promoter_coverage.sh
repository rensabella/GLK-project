## The Orthologous genes were found by blast, and the coverage behind genes were calcalated by bedtools multicov.

for s1 in AT NB SL OS ZM
do
for s2 in AT NB SL OS ZM
do
if [ "$s1" != "$s2" ]
then
bedtools multicov -bed <(awk '$4>0' promoter_site_${s1}_${s2}_blast|awk '$8>0'|cut -f3-5) -bams rmdup/${s1}_GLK1.bam rmdup/${s1}_GLK2.bam > promoter_coverage_${s1}_${s2}_${s1}
bedtools multicov -bed <(awk '$4>0' promoter_site_${s1}_${s2}_blast|awk '$8>0'|cut -f7-9) -bams rmdup/${s2}_GLK1.bam rmdup/${s2}_GLK2.bam > promoter_coverage_${s1}_${s2}_${s2}
fi
done
done
