This folder is for replicating the simulation results.

"loop_scenario1_and_2.txt" is for submitting batch jobs to HPC environment to get the results for Scenarios 1-2.
The boxplots can be obtained by running "scenario1_and_2_plot.R" after all the results are returned.

"loop_scenario3.txt" is for submitting batch jobs for Scenario 3. Running "scenario3_plot.R" can generate the corresponding boxplots.

External R packages in need include:
- MCMCpack
- glmnet
- doParallel
- foreach
- monomvn
- dimRed
- coda
- fastICA
- NMF
- igraph
- loe
- RSpectra
- RANN
- pROC
- ggsci
- gridExtra
- matrixcalc
- ggplot2
- reshape2
- grid