setwd("~/Desktop/research information/phylopredict/output")
Results <- readRDS("CrossValidationResults.RDS")
source("./Dataprep.R")
source("./ContMap2.R")
library('phytools')
install.packages("viridis")

#port over of required values of J and K from CrossValidationStudy.R
# Filters out NAs for all traits
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


test <- rbind(Results[[1]][[1]], Results[[1]][[2]])
test2 <- cbind(test, "trait"=traits[1])

#combining of the standard deviation, real value and p values
table <- list()
for(i in 1:length(traits)){
    .tmp <- do.call(rbind, Results[[i]])
    .tmp <- as.data.frame(apply(.tmp, 2, unlist)) #apply applies to the dataframe the function written over the given margin)
    table[[i]] <- cbind("trait"= traits[i], .tmp)
    #table[[i]] <- cbind(table[i], k$dat[1])
    table[[i]] <- cbind(table[[i]], "realvalue"=log(K[[i]][[1]]+0.01))
    table[[i]] <- cbind(table[[i]], "Stnd.Dev."=K[[i]][[2]])
    table[[i]] <- cbind(table[[i]], "p"=pnorm(table[[i]][,4],table[[i]][,2],table[[i]][,3] ))
}
#cbind(table, "trait"=traits[1]), function to cbind table data.
print(table)

par(ask=TRUE, mfrow=c(2,3)) #splits graphs into 2 rows by 3 columns
for(i in 1:length(table)){
  plot(table[[i]]$mu, table[[i]]$realvalue, main=table[[i]]$trait[1])
  abline(0,1)
  hist(table[[i]]$p)
} 

#function gives scatterplot and histogram for each of 86 traits
table2 <- list()
tmp2 <- list()
for(i in 1:86)
  {
    a <- table[[i]][[2]]
    b <- table[[i]][[4]]
    tmp2 <- ks.test(a,b)
    table2[[i]] <- do.call(data.frame, list("trait"=traits[i], "D"=tmp2$statistic, "pvalue"=tmp2$p.value))
}

#function for KS test
#do.call(rbind, table2)
dt <- do.call(rbind, table2)
#dtsorted <- dt[order(3)] #order(dt, na.last = TRUE, decreasing = FALSE, method = c("radix"))
#print(dt)  
#sort table so D values are from smallest to highest top down
#o <- order(dt$D)
data.table <- do.call(rbind, table2)
dtsorted <- arrange(dt, D, )
#function used for sorting values

#par(ask=TRUE, mfrow=c(2,3)) #splits graphs into 2 rows by 3 columns
#GOOD ONES - models worked well for these guys
plot(table[[which(traits=="GLUCOSE")]][,c(2,4)], main="Glucose", xlim=c(2,9), ylim=c(2,7), ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="GLUCOSE")]][,c(2,4)], labels=rownames(table[[which(traits=="GLUCOSE")]]))

plot(table[[which(traits=="ASPARTATEAMINOTRANSFERASE")]][,c(2,4)], main="Aspartateaminotransferase", ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="ASPARTATEAMINOTRANSFERASE")]][,c(2,4)], labels=rownames(table[[which(traits=="ASPARTATEAMINOTRANSFERASE")]]))

plot(table[[which(traits=="CHOLESTEROL")]][,c(2,4)], xlim=c(2,7), ylim=c(2,8),main="Cholesterol", ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="CHOLESTEROL")]][,c(2,4)], labels=rownames(table[[which(traits=="CHOLESTEROL")]]))

#BAD ONES - model performed poorly
plot(table[[which(traits=="PLATELETCOUNT")]][,c(2,4)], main="Platelet Count", ylab="Real Value", xlab="Predicted Value")
plot(table[[which(traits=="SODIUM")]][,c(2,4)], main="Platelet Count", ylab="Real Value", xlab="Predicted Value")

plot(table[[which(traits=="CORTISOL")]][,c(2,4)], main="Cortisol", ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="CORTISOL")]][,c(2,4)], labels=rownames(table[[which(traits=="CORTISOL")]]))

plot(table[[which(traits=="MAGNESIUM")]][,c(2,4)])
plot(table[[which(traits=="LIPASE")]][,c(2,4)], main="Lipase", ylab="Real Value", xlab="Predicted Value")

progesterone
plot(table[[which(traits=="AZUROPHILS")]][,c(2,4)], main="AZUROPHILS", ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="AZUROPHILS")]][,c(2,4)], labels=rownames(table[[which(traits=="AZUROPHILS")]]))

#MEH ONES - models had acceptable levels of prediction
plot(table[[which(traits=="CALCIUM")]][,c(2,4)], main="Calcium", ylab="Real Value", xlab="Predicted Value")

plot(table[[which(traits=="LYMPHOCYTES")]][,c(2,4)], main="Lymphocytes", xlim=c(-1,3), ylim=c(-2,4), ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="LYMPHOCYTES")]][,c(2,4)], labels=rownames(table[[which(traits=="LYMPHOCYTES")]]))

plot(table[[which(traits=="URICACID")]][,c(2,4)], main="Uricacid", ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="URICACID")]][,c(2,4)], labels=rownames(table[[which(traits=="URICACID")]]))

plot(table[[which(traits=="PHOSPHORUS")]][,c(2,4)], main="Phosphorous", ylab="Real Value", xlab="Predicted Value")
plot(table[[which(traits=="POTASSIUM")]][,c(2,4)], main="potassium", ylab="Real Value", xlab="Predicted Value")


plot(table[[which(traits=="IRON")]][,c(2,4)], main="Iron", xlim=c(0,8), ylab="Real Value", xlab="Predicted Value")
identify(table[[which(traits=="IRON")]][,c(2,4)], labels=rownames(table[[which(traits=="IRON")]]))
#need to fix this guy to move text over inside the graph boundary



#Set names on the difference between predicted and real value using phylogeny tip labels
#code for the contmaps
names(table) <- names(tds)
data2plot <- setNames(log(abs(table[['GLUCOSE']]$mu - table[['GLUCOSE']]$realvalue)/sqrt(table[['GLUCOSE']]$sigma),10), tds$GLUCOSE$phy$tip.label)
contMap2(tds$GLUCOSE$phy, data2plot, lims=c(-1, 0.5), ftype="off", leg.txt="Log10 Relative Error")

makeContMap <- function(trait, lims){
  data2plot <- setNames(log(abs(table[[trait]]$mu - table[[trait]]$realvalue)/sqrt(table[[trait]]$sigma),10), tds[[trait]]$phy$tip.label)
  contMap2(tds[[trait]]$phy, data2plot, lims=lims, ftype="off", leg.txt="Log10 Relative Error")
}
makeContMap("CORTISOL", lims=c(-1,0.2))





###outdated code for the contmaps
#data2plot <- setNames(table[[32]]$mu - table[[32]]$realvalue, tds$GLUCOSE$phy$tip.label)
contMap(tds$GLUCOSE$phy, data2plot, lims = c(-.4, .4), ftype="off") 
#data2plot2 <- setNames(table[[32]]$mu - table[[32]]$realvalue, tds$ASPARTATEAMINOTRANSFERASE$phy$tip.label)
contMap(tds$ASPARTATEAMINOTRANSFERASE$phy, data2plot, lims = c(-.4, .4), ftype="off") 
#data2plot3 <- setNames(table[[32]]$mu - table[[32]]$realvalue, tds$CHOLESTEROL$phy$tip.label)
contMap(tds$CHOLESTEROL$phy, data2plot, lims = c(-.4, .4), ftype="off") 

lim
range(data2plot)

##extra code, not really used for the graphs on the poster
par(mfrow=c(2,3))
lapply(traits[1:6], mynewfun)

lapply(traits[1:6], function(trait, ...){
  plot(table[[which(traits==trait)]][,c(2,4)], ...)
  #identify(table[[which(traits==trait)]][,c(2,4)], labels=rownames(table[[which(traits==trait)]]))
})

lapply(table, function(x) head(x,1))

mynewfun(table, "GLUCOSE")
#mynewfun(table, "GLUCOSE", xlim=c(-3,10), ylim=c(-2,10))
par(mfrow=c(2,6))
#traits[1:6]
#lapply(traits[1:6], mynewfun)
lapply(table, function(x) tail(x)+5) #
lapply(table, function(x) tail(x,1)+5)
lapply(table, function(x) tail(x,1))
lapply(table, function(x) head(x,1))
par(mfrow=c(2,3)) #sets number of rows and columns for graphs in plots
lapply(traits[1:6], mynewfun)

#tree, table, k, K,
#used to make the phylogenetic tree

#code used to allow for plotting of data points for each trait
mynewfun <- function(trait, ...){ #... means pass on arguments inside here to the next one
  plot(table[[which(traits==trait)]][,c(2,4)], ...)
  #identify(table[[which(traits==trait)]][,c(2,4)], labels=rownames(table[[which(traits==trait)]]))
}
#code that was used to look at multiple graphs at once
mynewfun <- function(table, trait){
  plot(table[[which(traits==trait)]][,c(2,4)])
  identify(table[[which(traits==trait)]][,c(2,4)], labels=rownames(table[[which(traits==trait)]]))
}

mynewfun <- function(trait, ...){ #... means pass on arguments inside here to the next one
  plot(table[[which(traits==trait)]][,c(2,4)], ...)
  #identify(table[[which(traits==trait)]][,c(2,4)], labels=rownames(table[[which(traits==trait)]]))
}
par(mfrow=c(2,3))
lapply(traits[1:6], mynewfun)
lapply(traits[1:6], function(trait, ...){
  plot(table[[which(traits==trait)]][,c(2,4)], ...)
  #identify(table[[which(traits==trait)]][,c(2,4)], labels=rownames(table[[which(traits==trait)]]))
})
lapply(table, function(x) head(x,1))


#viridis color palettes 
#wesanderson color palettes 
#phylopic - can use to search for species and find best matching taxon and species 
"want to describe the model we used"
#look at number of species for each trait
