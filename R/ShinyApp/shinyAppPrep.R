library(phytools)
library(shiny)
library(phylolm)
library(aRbor)
setwd("~/repos/phylopredict/R/ShinyApp")
dat <- read.csv("../../data/NORMALS.csv")
tree <- multi2di(read.tree("../../data/mamm.tre"))
tree$edge.length <- tree$edge.length/max(branching.times(tree))
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

dat$genspec <- unname(sapply(gsub(" ", "_", tolower(as.character(dat[[1]]))), simpleCap))
dat <- tbl_df(dat)


outliers <- function(x, iqrfactor){
  quants <- quantile(x)
  lowerq <- quants[2]
  upperq <- quants[4]
  iqr = upperq-lowerq
  x > iqr*iqrfactor+upperq | x < lowerq - iqr*iqrfactor

}
