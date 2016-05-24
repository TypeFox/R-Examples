# Function to calculate the division coefficients per group of birds and the recovery probability per area
# Author: Fraenzi Korner-Nievergelt, Swiss Ornithological Institute, www.vogelwarte.ch
# August 2009, R 2.9.1
###################################################################################################
dc <- function(N, recmatrix, group.names=NA, area.names=NA, start=NA){

# N: vector of the number of ringed birds per group
# recmatrix: matrix containing the number of re-encountered birds per group and area
# the rows of the matrix represent the bird groups, the columns represent the destination areas
###################################################################################################

if(nrow(recmatrix)==ncol(recmatrix)){
  nx <- recmatrix
  x <- solve(nx, N)
}

if(nrow(recmatrix)<ncol(recmatrix)) print("Number of groups must be at least the number of areas.")

if(nrow(recmatrix)>ncol(recmatrix)){
  ss <- function(x){
    x.mat <- matrix(x, nrow=nrow(recmatrix), ncol=length(x), byrow=TRUE)
    sums <- recmatrix*x.mat
    sums1 <- apply(sums, 1, sum)
    sum((sums1-N)^2)
  }
  ifelse(length(start)<2, init <- N[1]/recmatrix[1,]/2, init <- start)
  x <- optim(par=init, ss)$par
}

if(length(area.names)<2) area.names <- paste(rep("Area", ncol(recmatrix)), 1:ncol(recmatrix))
rec.probs <- matrix(1/x, nrow=1)
colnames(rec.probs) <- area.names

x.mat <- matrix(x, nrow=nrow(recmatrix), ncol=length(x), byrow=TRUE)
N.mat <- matrix(N, ncol=ncol(recmatrix), nrow=length(N))
mig.rates.obs <- x.mat*recmatrix/N.mat      # migration rates towards A for both groups
if(length(group.names)<2) group.names <- paste(rep("Group", nrow(recmatrix)), 1:nrow(recmatrix))
colnames(mig.rates.obs) <- area.names
rownames(mig.rates.obs) <- group.names

result <- list(rec.probs=rec.probs, division.coef=mig.rates.obs)
result
}
#### end of function################################################################################


