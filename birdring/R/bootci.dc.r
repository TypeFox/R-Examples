# Function to calculate the bootstrap confidence interval of the division coefficients and the reencounter probabilities
# Author: Fraenzi Korner-Nievergelt, Swiss Ornithological Institute, www.vogelwarte.ch
# August 2009, R 2.9.1
###################################################################################################

bootci.dc <- function(N, recmatrix, interval=0.95, R=1000, group.names=NA, area.names=NA){
# this function works only if the function dc() is loaded
# N: vector of the number of ringed birds per group
# recmatrix: matrix containing the number of re-encountered birds per group and area
# the rows of the matrix represent the bird groups, the columns represent the destination areas
###################################################################################################


if(length(group.names)<2) group.names <- 1:length(N)
if(length(area.names)<2) area.names <- 1:ncol(recmatrix)

nfound <- apply(recmatrix, 1, sum)
tmat <- cbind(N-nfound, recmatrix)
dat <- data.frame(group=rep(1:length(N), N),
  rec=rep(rep(c(0, area.names),length(N)), as.numeric(t(tmat))))

mig.rates <- matrix(nrow=R, ncol=length(area.names)*length(group.names))
gr <- rep(group.names, each=length(area.names))
are <- rep(area.names, length(group.names))
colnames(mig.rates) <- paste("G", gr, "A", are)

rec.probs <- matrix(nrow=R, ncol=length(area.names))
colnames(rec.probs) <- area.names
boot.rec.probs <- matrix(ncol=length(area.names), nrow=R)
boot.div.coef <- array(dim=c(length(group.names), length(area.names), R))

for(i in 1:R){
  dat.boot <- dat[sample(1:length(dat$group), replace=TRUE),]
  N.boot <- table(dat.boot$group)
  recmat.boot <- table(dat.boot$group, dat.boot$rec)[,-1][,area.names]
  dc.boot <- dc(N.boot, recmat.boot)
  boot.rec.probs[i,] <- dc.boot$rec.probs
  boot.div.coef[,,i] <- dc.boot$division.coef
}

lower <- function(x) quantile(x, probs=(1-interval)/2)
upper <- function(x) quantile(x, probs=interval+(1-interval)/2)

rec.probs.lower <- apply(boot.rec.probs, 2, lower)
rec.probs.upper <- apply(boot.rec.probs, 2, upper)
div.coef.lower <- apply(boot.div.coef, c(1,2), lower)
div.coef.upper <- apply(boot.div.coef, c(1,2), upper)
names(rec.probs.lower) <- area.names
names(rec.probs.upper) <- area.names
colnames(div.coef.lower) <- area.names
colnames(div.coef.upper) <- area.names
rownames(div.coef.lower) <- group.names
rownames(div.coef.upper) <- group.names

rec.probs.lower[rec.probs.lower<0]<-0
rec.probs.upper[rec.probs.upper>1]<-1

div.coef.lower[div.coef.lower<0]<-0
div.coef.upper[div.coef.upper>1]<-1

result.boot <- list(rec.probs.lower=rec.probs.lower,
  rec.probs.upper=rec.probs.upper,
  div.coef.lower=div.coef.lower,
  div.coef.upper=div.coef.upper)
result.boot
}
### end of function#############################################################

