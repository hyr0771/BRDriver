# BRDriver
BRDriver
(1) This package includes the Matlab scripts of BRDriver and five cancer datasets including BRCA, COAD, LUAD, PRAD and THCA.

(2) The user can start BRDriver with main_BRDriver.m. You can choose one of these five cancer datasets and specify the selected input dataset by setting the parameter “Data” in main_BRDriver.m. For example, we choose five patients’ tumor coding gene expression data and lncRNA expression data（only including 5000 genes）in BRCA as an example case and produce an example dataset 'exam'. The user can specify Data = 'exam' in main_BRDriver.m and run main_BRDriver.m to get a demo.

(3) The output are two resulting matrix of the predicted driver genes, one is the result of coding genes, and the other is the result of lncRNAs. In the resulting matrix, the first row is the name of patient samples and the genes in each column are the predicted driver genes for the patient sample.

(4) BRDriver was written and tested on Matlab with version "R2019a".
