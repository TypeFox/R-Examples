posteriorN <-
function(p, nf=0, maxN=1000, ci.int=0.95, plot=TRUE, dist=FALSE){
# p = probability that a killed animal is detected by a seacher
# nf = number of carcasses found
# maxN = maximal possible number of fatalities
# ci.int = size of the credible interval that should be calculated
# plot: posterior distribution is plotted if TRUE 
# dist: posterior distribution is given if  TRUE
#---------------------------------------------------------------
# track of changes:
# fk, 16.1. 2015: changed to use the function dbinom instead of choose and math
#---------------------------------------------------------------

N <- 0:maxN
 # version 1.0 - 1.3:
 # if(nf==0) pN1 <- p*(1-p)^(N-nf) fk: think this is unscaled!
 # denom <- sum(choose(N, nf) * (1-p)^(N-nf))
 # pN <- choose(N, nf)*(1-p)^(N-nf)/denom
 # pN <- c(rep(0, nf), pN)

 # from version 1.4 onwards 
  pCgN <-  dbinom(nf, size=N, prob=p)
  pN <- pCgN/sum(pCgN)

if(plot) plot(N, pN, type="h", lwd=5, lend="butt", xlab="Number of fatalities", ylab="Posterior density")
index <- cumsum(pN)<ci.int
indexLower <- cumsum(pN)<(1-ci.int)/2
indexUpper <- cumsum(pN)<1-(1-ci.int)/2
if(nf==0) interval <- c(nf, min(N[!index]))   
if(nf>0)  interval <- c(min(N[!indexLower]), min(N[!indexUpper])) 
if(interval[2]>maxN) cat("Upper limit of CI larger than maxN! -> increase maxN\n")
expected.median <- min(N[!cumsum(pN)<0.5])
#expected.mean <- sum(pN*N)
results <- list(interval=interval, expected=expected.median, HT.estimate=nf/p)
if(dist==TRUE) results <- list(interval=interval,  expected=expected.median, 
                               HT.estimate=nf/p, pN=pN)
results
}

# Wrapper, since in the version 1.0 of carcass the name of the function was posterior.N
posterior.N <- function(p, nf=0, maxN=1000, ci.int=0.95, plot=TRUE, dist=FALSE) {
  posteriorN(p=p, nf=nf, maxN=maxN, ci.int=ci.int, plot=plot, dist=dist)
}