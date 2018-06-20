setwd("~/Desktop/research information/phylopredict/R")
install.packages("treeplyr")
library(treeplyr)
# where we keep all our scripts, dont touch Raw Data, if changed, set as new data file)
# wd is working directory, getwd(), setwd()

sp360 <- readRDS("../data/species360.rds")

sp360$SPECIES
tree <- read.tree("../output/tetrapods.tre")
td <- make.treedata(tree,sp360)
summary(td)

#dplyer stuff
tmp <- filter(td, !is.na(CALCIUM), CALCIUM.N > 6) %>% select(., starts_with("CALCIUM")) %>% mutate(.,CALCIUM.log = log(CALCIUM))#allows us to use column names as if they were already saved variables
# is.na is to filter for data that is true for the na (N/A)
# 
