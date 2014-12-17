## # Building a mammal tree
## Load tree and data into R
#setwd("~/repos/phylopredict/R/")
setwd("./MammalTree")
require(ape)
require(aRbor)
require(devtools)
require(bayou)
#setwd("~/repos/bmr/R/MammalTree")
require(devtools)
#install_github("fmichonneau/rotl")
require(rncl)
require(rotl)
require(rjson)
require(geiger)
require(foreach)
require(doParallel)
source("./otTimetreeFunctions.R")

#taxa <- "mammals"
#td <- readRDS(paste("../../output/data/", taxa, ".rds", sep=""))
#td <- make.treedata(td$phy, td$dat, name_column=1)
#tree <- td$phy
#dat <- td$dat
#taxalist <- rownames(dat)
tree <- read.tree("../../data/mamm.tre")
taxalist <- tree$tip.label

## # Querying OpenTree
registerDoParallel(cores=12)
options("scipen"=100, "digits"=4)
tax <- foreach(i=1:length(taxalist)) %dopar% try(rotl::tnrs_match_names(taxalist[i], context="Mammals"))
tax <- do.call(rbind, tax)
genspec <- unname(sapply(tax[,2], function(x) paste(strsplit(x, split=" ")[[1]][1:2],collapse=" ")))
tax_unique <- tax[!duplicated(genspec),]
ottids <- as.character(tax_unique$ott_id) 

#ottids <- unique(ottids)


synth_tree <- rotl::tol_induced_subtree(ott_ids=ottids)
#pdf("mammalsynthtree.pdf", height=15)
#plot(synth_tree, cex=0.25)
#plot(tree, cex=0.25)
#dev.off()

## Find all source trees with our focal taxa in them, filter out trees without branch lengths. 
studyTrees <- findsourcetrees(ottids)      
possibleTrees <- getStudyTrees(studyTrees)
possibleTrees$summary

possibleTrees$summary[8,]

timetrees <- possibleTrees$result[which(possibleTrees$summary$timeUnits=="Myr")]
timetrees <- lapply(timetrees, function(x) x[[1]])
timetrees_treeId <- rownames(possibleTrees$summary[which(possibleTrees$summary$timeUnits=="Myr"),])
timetrees_studyId <- unlist(studyTrees[match(timetrees_treeId, studyTrees[,2]),1])
studyMeta <- lapply(timetrees_studyId, function(x) rotl::get_study_meta(study_id = x))
citations <- lapply(studyMeta, function(x) x$nexml[[2]])
print(citations)


problemtree <- list(studyId="pg_2688", treeId="tree6240")
ptree <- rotl::get_study_tree(problemtree$studyId, problemtree$treeId, text_format="newick", file="probtree.tre")
tree_string <- ptree
tree_string <- gsub( "\\[.*?\\]", "", tree_string)
tree_string <- gsub(" ", "_", tree_string)
tree_string <- gsub("'", "", tree_string, fixed=TRUE)
ptree <- read.newick(text=tree_string)
plot(ptree)
timetrees[[1]] <- ptree

cals <- synthCalibrations(synth_tree, timetrees)
ncals <- lapply(cals, nrow)
cals <- do.call(rbind, cals)
BE_cals <- synthCalibrations(synth_tree, list(tree))[[1]]
calid <- sapply(1:nrow(cals), function(x){ paste(sort(cals[x,4:5]),collapse="")})

## Look for duplicate calibrations and span the max and min of that calibration for the age. 
calibrations <- NULL
for(i in 1:length(unique(calid))){
  tmp <- subset(cals, calid==unique(calid)[i])
  if(nrow(tmp)>1){
    minage <- min(tmp[,2])
    maxage <- max(tmp[,3])
    tmp[1,2:3] <- c(maxage, minage)
    rowi <- tmp[1,]
  } else {
    rowi <- tmp
  }
  calibrations <- rbind(calibrations, rowi)
}

## Pathd8 first with BE
target <- synth_tree
target$tip.label <- unname(sapply(synth_tree$tip.label,function(x) strsplit(x, "_ott")[[1]][1]))
target$edge.length <- rep(1,nrow(target$edge))
target <- chronopl(target, lambda=1)
tr <- replaceHashes(target, BE_cals)
phypd8 <- PATHd8.phylo(tr$scion, tr$table, base = ".tmp_PATHd8", rm = FALSE)

##Pathd8 again with better cals:
#phypd8$node.label <- target$node.label
## A set of hand-entered calibrations from Meredith 2011 supplement:
#meredithCalibrations <- list(
#                              data.frame("MRCA"=NA, "MaxAge"=238.16, "MinAge"=210.71, "taxonA"="Ornithorhynchus_anatinus", "taxonB"="Hipposideros_galeritus"),#All mammals
#                              data.frame("MRCA"=NA, "MaxAge"=116.81, "MinAge"=92.12, "taxonA"="Cercopithecus_mitis", "taxonB"="Cyclopes_didactylus"),#Placentalia
#                              data.frame("MRCA"=NA, "MaxAge"=82.64, "MinAge"=79.07, "taxonA"="Elephantulus_rufescens", "taxonB"="Procavia_capensis"),#Afrotheria
#                              data.frame("MRCA"=NA, "MaxAge"=68.42, "MinAge"=62.66, "taxonA"="Chaetophractus_vellerosus", "taxonB"="Cyclopes_didactylus"),#Xenarthra
#                              data.frame("MRCA"=NA, "MaxAge"=67.39, "MinAge"=65.22, "taxonA"="Paranyctimene_raptor", "taxonB"="Molossus_molossus"),#Chiroptera
#                              data.frame("MRCA"=NA, "MaxAge"=82.78, "MinAge"=77.46, "taxonA"="Panthera_leo", "taxonB"="Manis_javanica"),#Carnivora+Philodota
#                              data.frame("MRCA"=NA, "MaxAge"=72.95, "MinAge"=66.03, "taxonA"="Thomomys_bottae", "taxonB"="Glis_glis"),#Rodentia
#                              data.frame("MRCA"=NA, "MaxAge"=87.00, "MinAge"=75.40, "taxonA"="Cercopithecus_mitis", "taxonB"="Glis_glis"),#Glires
#                              data.frame("MRCA"=NA, "MaxAge"=90.49, "MinAge"=78.11, "taxonA"="Cercopithecus_mitis", "taxonB"="Galago_senegalensis"),#Primates+Dermoptera
#                              data.frame("MRCA"=NA, "MaxAge"=90.27, "MinAge"=71.68, "taxonA"="Didelphis_virginiana", "taxonB"="Macrotis_lagotis"),#Marsupialia
##                              data.frame("MRCA"=NA, "MaxAge"=49.57, "MinAge"=30.36, "taxonA"="Ornithorhynchus_anatinus", "taxonB"="Tachyglossus_aculeatus")#Monotremata
#                              )
meredith <- read.tree("../../data/meredith-mammals.tre")
meredithotts <- rotl::tnrs_match_names(meredith$tip.label, context_name = "Mammals")$ott_id
meredithtax <- gbresolve(phypd8$tip.label)

meredithCalibrations <- congruify.phylo(meredith, phypd8, meredithtax, tol=0, scale=NA)$calibrations

calibrations <- rbind(calibrations, meredithCalibrations)
tr <- replaceHashes(phypd8, calibrations)
#excluded <- c(144,146,147,148,149, 151,152,153,154,171)
rcalibrations <- reconcileCalibrations(tr$scion, tr$table)

finalTree <- PATHd8.phylo(tr$scion, rcalibrations, base = ".tmp_PATHd8", rm = FALSE)
finalTree <- reorder(finalTree, "postorder")

BETree <- reorder(tr$scion, "postorder")

par(mfrow=c(1,2))
plot(finalTree, show.tip.label=FALSE)
abline(v=max(nodeHeights(finalTree))-65, lty=2)
plot(BETree, show.tip.label=FALSE)
abline(v=max(nodeHeights(BETree))-65, lty=2)

write.tree(finalTree, "../../output/fullnewmammal.tre")
