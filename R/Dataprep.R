#rm(list=ls(all=TRUE)) #erases everything from the environment 
setwd("~/Desktop/research information/phylopredict/R")
#install.packages("treeplyr")
library(treeplyr) #loaded treeplyr, if function is not found, package not loaded, which is what library is for library(treeplyr) -> better for script, require(treeplyr) -> does not give error if treeplyr not found) 
# where we keep all our scripts, dont touch Raw Data, if changed, set as new data file)
# wd is working directory, getwd(), setwd()
sp360 <- readRDS("../data/species360.rds")
tree <- read.tree("../output/tetrapods.tre")
td <- make.treedata(tree,sp360)
summary(td)

#dplyer stuff
traits <- unique(colnames(td$dat))
traits <- traits[1:86]
tds <- list() #list is group of unrelated stuff; important

for(i in traits){
  .td <- select(td, starts_with(i))
  .td <- filter(.td, !is.na(.td$dat[[1]]), .td$dat[[1]] != 0)
  .td <- mutate(.td, trait.log = log(.td$dat[[1]]))
  tds[[i]] <- .td
  print(tds[i])
}
?starts_with
#tmp <- filter(td, !is.na(CALCIUM)) %>% select(., starts_with("CALCIUM")) %>% mutate(.,CALCIUM.log = log(CALCIUM))#allows us to use column names as if they were already saved variables
# is.na is to filter for data that is true for the na (N/A)
# 
#sp360$SPECIES
#tree <- read.tree("../output/tetrapods.tre")
#td <- make.treedata(tree,sp360)
#summary(td)

#tree <- NULL
#read.tree("../output/tetrapods.tre")

#sp360 <- readRDS("../data/species360.rds")

#make for loop for each trait where species had missing values, then make table where each trait is split into its own tree species was matched with 