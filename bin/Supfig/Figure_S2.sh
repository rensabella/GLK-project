## Calculate and plot command (Arabidopsis)                             
# bigwig2bedGraph                               
ls *bw | while read id; do label=$(basename $id ".bw"); ./bigWigToWig ${id} ${label}.bedGraph; done            

# bedGraph/FRiP (python script)                         
python bedGraph_FRiP.py 28.41 AtGLK1_rep1.bedGraph AtGLK1_rep1_FRiP.bedGraph &
python bedGraph_FRiP.py 26.99 AtGLK1_rep2.bedGraph AtGLK1_rep2_FRiP.bedGraph &
python bedGraph_FRiP.py 33.25 AtGLK2_rep1.bedGraph AtGLK2_rep1_FRiP.bedGraph &
python bedGraph_FRiP.py 32.88 AtGLK2_rep2.bedGraph AtGLK2_rep2_FRiP.bedGraph &

# sort bedGraph                         
ls *_FRiP.bedGraph|while read id; do label=$(basename $id ".bedGraph"); sort -k1,1 -k2,2n ${id} > ${label}.sortBedGraph; done

# bedGraph2bigwig                               
ls *.sortBedGraph | while read id; do label=$(basename $id ".sortBedGraph"); ./bedGraphToBigWig ${id}_FRiP.sortBedGraph TAIR10.chromsize ${id}_FRiP.bigwig; done

# calculate matrix                              
computeMatrix scale-regions -S \
AtGLK1_rep1_FRiP.bigwig \
AtGLK1_rep2_FRiP.bigwig \
AtGLK2_rep1_FRiP.bigwig \
AtGLK2_rep2_FRiP.bigwig \
-R At_G1G2_overlap.bed At_G1_only.bed At_G2_only.bed \
--regionBodyLength 150 \
--beforeRegionStartLength 1000 \
--afterRegionStartLength 1000 \
--binSize 10 \
--numberOfProcessors 10 \
-o At_GLK.profile.mat.gz \
--samplesLabel AtGLK1_1 AtGLK1_2 AtGLK2_1 AtGLK2_2

# plot                          
plotProfile -m At_GLK.profile.mat.gz \
-out At_GLK.profile.pdf \
--plotTitle At_GLK --perGroup \
--startLabel start \
--endLabel end
