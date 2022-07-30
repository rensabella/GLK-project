# Find motifs in top 2000 peaks of 5 species by HOMER 
cd ../../data/motif_find_fasta
for s in AT SL NB OS ZM
do
findMotifs.pl ${s}GLK1.fasta fasta HOMER_${s}_GLK1 -p 40
findMotifs.pl ${s}GLK2.fasta fasta HOMER_${s}_GLK2 -p 40
done

