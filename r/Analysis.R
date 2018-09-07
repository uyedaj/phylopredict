source("./Dataprep.R")
library("phylolm")


rownames(tds[[1]]$dat) <- tds[[1]]$phy$tip.label
res <- phylolm(trait.log~1, data = tds[[1]]$dat, phy = tds[[1]]$phy, model = "BM") #formula, variable saying that y = mx+b, single value sets y = b)
#sigma 2 is rate of evolution, size of drunkards walks
# intercept is the b, starting value of the branches, the root value, lots of uncertainty

#put this in a loop so that you have a list of all the traits and save it
#fit models to every trait and loop it over save output as list 
 #.td <- select(td, starts_with(i))
    #.td <- filter(.td, !is.na(.td$dat[[1]]), .td$dat[[1]] != 0)
    #.td <- mutate(.td, trait.log = log(.td$dat[[1]]))
    #tds[[i]] <- .td



  #print(inter[i])
  #print(sig2[[i]]$sigma2)
  #print(ads[[i]]$sigma2)
#cbind(traits, sig2, inter)
#data.table <- cbind(traits, inter, sig2)
#ads[[1]]$sigma2

### clean code begins here
# data for Brownian Motion
adsBM <- list()
sig2BM <- list()
#interBM <- rep(NA, length(traits)) also works in place of list
interBM <- list()
alphaBM <- list()
rBM <- list()
for(i in 1:length(traits)){
  res <- phylolm(trait.log~1, data = tds[[i]]$dat, phy = tds[[i]]$phy, model = "BM") #formula, variable saying that y = mx+b, single value sets y = b)
  adsBM[[i]] = res
  sig2BM[[i]] <- adsBM[[i]]$sigma2
  interBM[i] <- adsBM[[i]]$coefficients
  alphaBM[i] <- NA
  rBM[i] <- NA
}
#data.table <- cbind(traits, sig2BM, interBM, alphaBM, rBM)
#print(data.table)

# data for OU
adsOU <- list()
sig2OU <- list()
interOU <- list()
alphaOU <- list()
rOU <- list()
for(i in 1:length(traits)) {
  res <- phylolm(trait.log~1, data = tds[[i]]$dat, phy = tds[[i]]$phy, model = "OUrandomRoot")
  adsOU[[i]] = res
  sig2OU[[i]] <- adsOU[[i]]$sigma2
  interOU[i] <- adsOU[[i]]$coefficients
  alphaOU[i] <- adsOU[[i]]$optpar
  rOU[i] <- NA # use this for table, numeric type for N/A, which is character
}
#data for EB 
adsEB <- list()
sig2EB <- list()
interEB <- list()
alphaEB <- list()
rEB <- list()
for(i in 1:length(traits)) {
  res <- phylolm(trait.log~1, data = tds[[i]]$dat, phy = tds[[i]]$phy, model = "EB")
  adsEB[[i]] <- res
  sig2EB[[i]] <- adsEB[[i]]$sigma2
  interEB[i] <- adsEB[[i]]$coefficients
  rEB[i] <- adsEB[[i]]$optpar
  #rOU[i] <- adsEB[[i]]$optpar #r is stored in optpar of OU parameter
}

#combined data table of BM, OU and EB)

data.table <- cbind(traits, sig2BM, sig2OU, sig2EB, interBM, interOU, interEB, alphaBM, alphaOU, alphaEB, rBM, rOU, rEB)
print(data.table)

# fix "N/A" to NA 
# do rbind, make 3 separate tables, with same column names, then rbind 


  
