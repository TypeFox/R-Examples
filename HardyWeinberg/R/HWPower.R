HWPower <- function(n=100,nA=100,pA=0.5,y=c(AA=25,AB=50,BB=25),alpha=0.05,theta=4,f=NULL,test="exact",alternative="two.sided", pvaluetype="selome",cc=0.5) {
if(!missing(pA)) {
    if((pA > 1) | (pA < 0)) stop("Allele frequency not in the range (0,1)")
    nA <- round(2*n*pA,digits=0)
}
if(!missing(y)) {
    n <- sum(y)
    nA <- mac(y)
}
if(!is.wholenumber(n) | !is.wholenumber(nA)) {
  warning("n or nA are not integers, and will be rounded")
  n <- round(n,digits=0)
  nA <- round(nA,digits=0)
}
if((n <=0) | (nA < 0)) stop("n or nA are non-positive")
nB <- 2*n-nA
if(nB < nA) {
  warning("nA is not the minor allele count, using 2n-nA instead")
  nA <- nB
  nB <- 2*n-nA
}
if(!is.null(f)) {
   pA <- nA/(2*n)
   pB <- 1-pA
   minaf <- min(pA,pB)
   fmin <- -minaf/(1-minaf)
   if(f < fmin | f > 1) stop("HWPower: inbreeding coefficient out of range")
   theta <- FtoTheta(pA,f)
}
if(nA%%2 == 1) nAB <- seq(1,nA,by=2)
if(nA%%2 == 0) nAB <- seq(0,nA,by=2)
nAA <- (nA-nAB)/2
nBB <- (nB-nAB)/2
X <- cbind(nAA,nAB,nBB)
colnames(X) <- c("AA","AB","BB")
probs <- HWExactPower(X[1,],theta=theta)$prob
p.value.vec <- NULL
for (i in 1:nrow(X)) {
   p.value <- switch(test,
                     exact = HWExactPower(X[i,],alternative=alternative,pvaluetype=pvaluetype,verbose=FALSE,theta=4)$pval,
                     chisq = HWChisq(X[i,],cc=cc,alpha)$pval,
                     stop("invalid value for parameter test"))
   p.value.vec <- c(p.value.vec,p.value)
}
ind.sig <- p.value.vec < alpha
pw <- sum(probs[ind.sig])
return(pw=pw)
}
