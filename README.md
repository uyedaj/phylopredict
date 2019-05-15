# phylopredict repository

This repository holds code for conducting a Cross Validation study of phylogenetic prediction from the species 360 database. 

The key scripts for generating this data are in the ~/R/ directory and include: 

R/Dataprep.R (prepares data to be analyzed)

R/CrossValidationStudy.R (runs cross validation study across multiple cores)

R/CrossValidationTable.R (collates results and evaluates prediction performance)

output/fullCrossValidationResults.rds (Final results for all traits)
