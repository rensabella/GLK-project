# Cross-species comparison of transcription factor binding revealed cistrome plasticity during plant evolution

This is an introduction about the code we used in this project. The scripts are all in the folder [bin](https://github.com/rensabella/GLK-project/tree/main/bin "bin").


## 1. GLK binding motif prediction

At first we detected the novel motifs just with [HOMER](https://github.com/rensabella/GLK-project/blob/main/bin/GLK_binding_motif/HOMER_find_motif.sh). To overcome the limitation of this software, we developed an [K-mer model](https://github.com/rensabella/GLK-project/blob/main/bin/GLK_binding_motif/kmer_model.ipynb) to find motifs contribute most in the binding sites.  The input format could refer to the [tables](https://github.com/rensabella/GLK-project/tree/main/data/kmer_data) we uploaded. These sequences are extracted from the peak summit regions.

## 2. Species-specific analysis
The GLK target genes we identified could be divided into five groups according to their conservation in five species. This [script](https://github.com/rensabella/GLK-project/blob/main/bin/Species_specific/species_specific_analysis.R) is to find the conservation of a target.
## 3. Function analysis

Uploaded the GO result files from agriGO and plot the results. 

## 4. Duplication pattern

Find the pattern in maize genome duplication.
## 5. DAP model
First, the [comparison of models](https://github.com/rensabella/GLK-project/blob/main/bin/DAP_model/compare%20models.ipynb) with multiple features shows each of our selected features have important effects on the model. In the [best model script](https://github.com/rensabella/GLK-project/blob/main/bin/DAP_model/Best%20model%20plot.ipynb), we resample the training and test dataset for 500 times to avoid errors caused by small dataset, which achieves high accuracy. 
