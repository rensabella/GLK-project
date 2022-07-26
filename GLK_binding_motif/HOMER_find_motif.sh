# Find motifs in top 2000 peaks of 5 species by HOMER 
for s in AT SL NB OS ZM
do
cd $s
findMotifs.pl ${s}GLK1.fasta fasta HOMER_${s}_GLK1 -p 40
findMotifs.pl ${s}GLK2.fasta fasta HOMER_${s}_GLK2 -p 40
cd ..
done

