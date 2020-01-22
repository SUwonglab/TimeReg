# TimeReg
time course regulatory analysis from paired gene expression and chromatin accessibility time course data 

At each time point we use the expression and accessibility data to define two types of scores for regulatory relations. 
The trans-regulatory score (TRS) quantifies the regulatory strength of a transcription factor (TF) on a target gene (TG), 
while the cis-regulatory score (CRS) quantifies the regulatory strengths of a regulatory element (RE) on a target gene. 
Based on these scores, we use Non-negative matrix factorization to extract the core regulatory modules that characterize different 
biological processes and/or subpopulations of cells, and we identify diver TFs (i.e. TFs driving expression changes between adjacent
time points) based on changes in TFs expression and changes in CRS on REs with increasing accessibility.

#Install
https://github.com/SUwonglab/TimeReg/archive/master.zip
