# Kanari et al. Leukemia 2025
This repository contains code and links to raw data used for: *A three-dimensional ex vivo model recapitulates in vivo features and drug resistance phenotypes in childhood acute lymphoblastic leukemia*

# Dataset composition
1. The 3D data include leukemic, endothelial and mesenchymal stromal cells cultured in the Ectica plates available at:
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282806

2. The *in vivo* data include leukemic cells from NSG mice at full PDX engraftment available at:
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284242 

## Analysis 
### Preprocessing leukemia cells
The procesing and integration and downstream analysis of the leukemic PDX data can be found in the [leukemic_cells_paper.R](https://github.com/mkanari/scRNAseq_analysis/blob/main/leukemia_cells_paper.R)

### Preprocessing mesenchymal stromal and endothelial cells
The processing and downstream analysis of the mesenchymal stromal and endothelial cells can be found in the [supporing_cells_paper.R](https://github.com/mkanari/scRNAseq_analysis/blob/main/supporting_cells_paper.R)

### Figures
The graphs for the figures were created with the code found in [figures_paper.R](https://github.com/mkanari/scRNAseq_analysis/blob/b5ee7778aeaca6fc446d70967c4ce65cea1f90a6/figures_paper.R) 

### Supplementray figures
The graphs for the supplementary figures were created with the code found in [suppl_figures_paper.R
](https://github.com/mkanari/scRNAseq_analysis/blob/main/suppl_figures_paper.R)

### Session and packages information
The [session info](https://github.com/mkanari/scRNAseq_analysis/blob/main/session_info.md) and [packages information](https://github.com/mkanari/scRNAseq_analysis/blob/main/installed_packages.csv) can be found in the respective files. 

