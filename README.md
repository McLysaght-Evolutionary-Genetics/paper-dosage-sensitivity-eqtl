# Evolution of dosage-sensitive genes by tissue-restricted expression changes

## Alan M. Rice, Yuanshuo Li, Pauric Donnelly, and Aoife McLysaght

This repository contains code and data used for analyses in the manuscript.

### Code tested on

* Python - *v3.6.6*
* Bedtools - *v2.27.1*
* Bedops - *v2.4.12*
* R - *v3.5.1*

Python dependancies:
* pandas: 0.23.3
* numpy: 1.15.0
* scipy: 1.1.0
* mysql-connector-python: 2.0.4
* sqlalchemy: 1.2.10

R dependancies:
* forcats: 0.3.0
* GeneOverlap: 1.16.0
* ggplot2: 3.0.0
* IRdisplay: 0.5.0
* IRkernel: 0.8.12.9000
* RColorBrewer: 1.1-2
* reshape2: 1.4.3
* scales: 1.0.0
* viridis: 0.5.1


### Important files to note

* `analysis/`
  * `notebooks/`
    * `prepareData`
      * `eQTLImportAndCorrect.py3.ipynb` - Import GTEx v7 eQTL data into MySQL notebook and correct for multiple testing
      * `GTExExpressedAndTestedGenes.python.ipynb` - Identify genes expressed and tested for eQTLs in GTEx v7
      * `HumanOhnologsSSDsSingletons.py3.ipynb` - Identify human ohnologs, singletons, and SSDs
      * `HumanCNVRGenes.py3.ipynb` - Identify human genes in Zarrei et al.
      * `ExACpLIScore.py3.ipynb` - ExAC genes with high and low pLI scores
      * `eQTLsBroadTissueBreadth.py3.ipynb` - Identify broad tissue breadth eQTLs
      * `countSNPsTestedPerGene.py3.ipynb` - Count number of SNPs tested per gene

### Citation

To follow...
