testGoF <-
function(data, type, numSims=10, plot=TRUE,  main){
if(missing(data) || missing(type))
stop("data and/or type is missing.")

if(numSims <= 0)
stop("The number of simulations must be an integer greater than 0.")

numGraphs <- ncol(data)
gstar <- estGStar(data)
tau <- estTau(data, type, gstar)
nodes <- getNumNodes(data, type)
edges <- getNumEdges(nodes, type)

#Calculate the prob for a graph k-distance away and calculate the expected counts of trees k-distance from gstar
possEdges <- 0:edges
distProb <- exp(lgamma(edges+1) - lgamma(edges-possEdges+1) - lgamma(possEdges+1) - tau*possEdges-edges * log(1+exp(-tau)))
expCounts <- cbind(possEdges, numGraphs*distProb)

#Calculate the observed counts of trees k-distance from gstar
distToGStar <- apply(data, 2, function(x, g, type){calcDistance(x, g, type)}, g=as.vector(gstar), type=type)
distTable <- table(distToGStar)
obsCounts <- cbind(as.integer(names(distTable)), as.integer(distTable))

#Combine observed and expected counts into a single object
obsExpCounts <- merge(expCounts, obsCounts, by=1, all=TRUE)
colnames(obsExpCounts) <- c("dist", "expected", "observed")
obsExpCounts[is.na(obsExpCounts)] <- 0
methodA <- "Chisq"

#Combine distances with a theoretical count < 5 (1 if there aren't enough at 5)
signifExpect <- obsExpCounts$expected >= 5 
if(sum(signifExpect) <= 1){
signifExpect <- obsExpCounts$expected >= 1
if(sum(signifExpect) < 1){
stop("Expected counts below the threshold of 1. Goodness of fit cannot be calculated")
}else{
methodA <- "MC"
warning("Expected counts below the threshold of 5; threshold lowered to 1 \n 
P-value computed using Monte-Carlo simulation instead of asymptotic distribution.")
}
}

#Determine main concentration of graphs
sumSignifExp <- sum(signifExpect)
lowerBound <- which(cumsum(signifExpect)==1)
upperBound <- length(signifExpect)-which(cumsum(rev(signifExpect))==1)+1
signifExpect[c(lowerBound, upperBound)] <- FALSE 

#Collapse groups below/above the bounds
lowerCount <- colSums(obsExpCounts[1:lowerBound, c("expected", "observed")])
upperCount <- colSums(obsExpCounts[upperBound:length(signifExpect), c("expected", "observed")])
mainCount <- obsExpCounts[signifExpect, c("expected", "observed")]

#Bind all groups into 1 table
oTable <- rbind(lowerCount, mainCount, upperCount)
rownames(oTable) <- c(paste("<=", obsExpCounts[lowerBound, 1]), obsExpCounts[signifExpect, 1], paste(">=", obsExpCounts[upperBound, 1]))

#Compute Pearson Chi-Squared Statistics
pearsonStats <- sum(((oTable$observed-oTable$expected)^2)/oTable$expected)
df <- nrow(oTable)-1

#Compute the Chi-Squared Statistics
chisq <- chisq.test(oTable$observed, p=oTable$expected/numGraphs, simulate.p.value=ifelse(methodA=="MC",TRUE,FALSE), B=numSims)
pvalueC <- chisq$p.value

#Compute the G Statistics if needed
if(sum(oTable$observed==0)==0){
gStats <- 2*sum(oTable$observed * log(oTable$observed/oTable$expected))
gdf <- nrow(oTable)-1
if(gdf <= 0){
gdf <- nrow(oTable)
warning("df is zero or negative; df replaced by the number of cells. (Conservative Test)")
}
pvalueG <- pchisq(gStats, df=gdf, ncp=0, lower.tail=FALSE, log.p=FALSE)
}else{
gStats <- NA
pvalueG <- NA
}

pval <- ifelse(is.na(pvalueG), pvalueC, pvalueG)
ptype <- ifelse(is.na(pvalueG), ifelse(methodA=="MC", "Monte-Carlo simulation", "Pearson Chi-square"), "G-test Statistics")

#Plot the observed vs expected data
if(plot){
mycolor <- c("red", "blue")
mylegend <- c("Expected", "Observed")
signifExpect[c(lowerBound, upperBound)] <- TRUE 
dist <- obsExpCounts[signifExpect, 1]
matplot(dist, oTable, pch=19, type="p", col=mycolor)
legend("topright", legend=mylegend, col=mycolor, pch=19)
if(missing(main))
title(c(paste(ptype), paste("P-value:", round(pval, 2))))
else
title(main)
}

results<- list(ptype, df, pval, oTable)
names(results) <- c("Method", "df", "pvalue",  "table")
return(results)
}
