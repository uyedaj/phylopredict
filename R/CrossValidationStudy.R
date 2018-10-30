#environment setup
setwd("~/Desktop/research information/phylopredict/R")
source("./Dataprep.R")
library(treeplyr)
library(phylolm)

# Loading data
sp360 <- readRDS("../data/species360.rds")
data.table <- as.data.frame(readRDS("./coefficients.rds"))
tree <- read.tree("../output/tetrapods.tre")
td <- make.treedata(tree,sp360)
summary(td)

# Defining function 
predictMissingTips <- function (tree, dat, sig2, root) {
  V <- phytools::vcvPhylo(tree, anc.nodes = FALSE)
  X <- dat
  unknown <- is.na(X)
  known <- !unknown
  Vkk <- V[known, known]
  Vuu <- V[unknown, unknown]
  Vku <- V[known, unknown]
  Vuk <- V[unknown, known]
  iVkk <- solve(Vkk)
  sigmabar <- as.matrix(Matrix::forceSymmetric(Vuu - Vuk %*% 
                                                 iVkk %*% Vku))
  cholSigmabar <- chol(sigmabar)
  mubarmat <- Vuk %*% iVkk
  prevalues <- list(V = V, X = X, unknown = unknown, known = known, 
                    Vkk = Vkk, Vuu = Vuu, Vku = Vku, Vuk = Vuk, iVkk = iVkk, 
                    sigmabar = sigmabar, mubarmat = mubarmat, cholSigmabar = cholSigmabar)
  X <- prevalues$X
  Vuk <- sig2 * prevalues$Vuk
  iVkk <- (1/sig2) * prevalues$iVkk
  Vku <- sig2 * prevalues$Vku
  Vuu <- sig2 * prevalues$Vuu
  known <- prevalues$known
  unknown <- prevalues$unknown
  mu <- rep(root, length(X))
  muk <- mu[known]
  muu <- mu[unknown]
  mubar <- t(muu + Vuk %*% iVkk %*% (X[known] - muk))
  sigmabar <- Vuu - Vuk %*% iVkk %*% Vku
  #res <- MASS::mvrnorm(1, mubar, sigmabar)
  #pars.new <- pars
  #pars.new$missing.pred <- res
  #hr = Inf
  #type = "impute"
  return(list(mu=mubar, sigma=sigmabar))
}

# Estimating BM parameters for the model
rownames(tds[[1]]$dat) <- tds[[1]]$phy$tip.label
res <- phylolm(trait.log~1, data = tds[[1]]$dat, phy = tds[[1]]$phy, model = "BM")
  adsBM <- list()
  sig2BM <- list()
  interBM <- list()  #interBM <- rep(NA, length(traits)) also works in place of list
  alphaBM <- list()
  rBM <- list()
    for(i in 1:length(traits))
    {
      res <- phylolm(trait.log~1, data = tds[[i]]$dat, phy = tds[[i]]$phy, model = "BM") #formula, variable saying that y = mx+b, single value sets y = b)
      adsBM[[i]] = res
      sig2BM[[i]] <- adsBM[[i]]$sigma2
      interBM[i] <- adsBM[[i]]$coefficients
      alphaBM[i] <- NA
      rBM[i] <- NA
    }

# Filter out NAs for all traits
traits <- unique(colnames(td$dat))
traits <- traits[1:86]
  K <- list()
  J <- list()
    for(i in traits) 
    {
      #i = traits[1] #string has to match with the starts_with string. column names must match up
      k <- select(td, starts_with(i))
      K[[i]] <- filter(k, !is.na(k$dat[[1]]), k$dat[[1]] !=0)
      #print(K[[i]])
    }
  
#code for multicore use
library(foreach)
library(doParallel)  
registerDoParallel(cores = 50)

# Cross validation study    
missing <- list()
missing_trait <- list()
    for (i in 1:length(K)) 
      {
      missing_trait[[i]] <- list()
      fulldata <- log(K[[i]]$dat[[1]]+0.01)
      missing_trait[[i]] <- foreach(j=1:nrow(K[[i]]$dat)) %dopar%
        { 
          tmpdata <- fulldata
          tmpdata[j] <- NA 
          missing <- predictMissingTips(K[[i]]$phy, tmpdata, sig2=sig2BM[[i]], root=interBM[[i]])
          missing
          #missing_trait[[i]][[j]] <- missing
          #do.call(rbind, missing_trait)
          #data.table <- missing_trait
        }
      }
#how to save as a file  
saveRDS(missing_trait,"../output/CrossValidationResults.RDS")

 