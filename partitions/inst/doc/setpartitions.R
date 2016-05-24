### R code from vignette source 'setpartitions.Rnw'

###################################################
### code chunk number 1: load.library
###################################################



###################################################
### code chunk number 2: setpartitions.Rnw:82-83
###################################################
require(partitions)


###################################################
### code chunk number 3: alleles
###################################################
setparts(c(2,2,1))


###################################################
### code chunk number 4: use_achims_idea
###################################################
a <- c(1,2,3,2,1)
split(seq_along(a),a)


###################################################
### code chunk number 5: setpartitions.Rnw:362-365
###################################################
options(width=63)
m <- 9
n <- 3


###################################################
### code chunk number 6: setpartsrestrictedparts
###################################################
jj <- setparts(restrictedparts(m,n,include.zero=FALSE))
summary(jj)


###################################################
### code chunk number 7: definepenaltyfunctionF
###################################################
tau <- 1:9
slowest <-  function(x) max(tapply(tau, x, sum))


###################################################
### code chunk number 8: show_the_best
###################################################
time.taken <- apply(jj,2,slowest)
minimal.time <- sum(tau)/n
jj[,time.taken == minimal.time]


###################################################
### code chunk number 9: define.E.matrix
###################################################
E <- matrix(c(0,0,0,1,1,0,0,0,1,1,1,1,1,0,0), 5,3)
dimnames(E) <-  list(
                     evidence=paste("E",1:5,sep=""),
                     crime=paste("C",1:3,sep="")
                     )


###################################################
### code chunk number 10: print.E
###################################################
E


###################################################
### code chunk number 11: set.evidence.matrix
###################################################
a <- 
structure(
          c(2, 2, 4, 4, 3, 1, 2, 1, 1, 2, 2, 1, 4, 2, 1, 2, 1, 
            4, 2, 4, 2, 1, 4, 4, 1, 2, 1, 4, 2, 1, 1, 2, 1, 1, 2),
          .Dim = c(5L, 7L),
          .Dimnames = list(
            evidence=paste("E",1:5,sep=""),
            crime=paste("C",1:7,sep="")
            )
)


###################################################
### code chunk number 12: show.a
###################################################
a


###################################################
### code chunk number 13: calc.likelihoods
###################################################
genbeta <- function(x){gamma(sum(x))/prod(gamma(x))}

lp <- function(evidence, alpha, partition){
  r <- length(unique(partition))
  evidence <- as.factor(evidence)
  levels(evidence) <- seq_along(alpha)
  out <- 1
  for(k in unique(partition)){  #k is a perp
    thisperp <- partition==k
    evidence.thisperp <- evidence[thisperp]
    no.of.crimes.thisperp <- sum(thisperp)
    out <- out/genbeta(alpha+table(evidence.thisperp))
  }
  return(out*genbeta(alpha)^r)
}

sp <- setparts(7)
l1 <- apply(sp,2, function(x){lp(evidence=a[1,],alpha=rep(1,2),partition=x)})
l2 <- apply(sp,2, function(x){lp(evidence=a[2,],alpha=rep(1,2),partition=x)})
l3 <- apply(sp,2, function(x){lp(evidence=a[3,],alpha=rep(1,4),partition=x)})
l4 <- apply(sp,2, function(x){lp(evidence=a[4,],alpha=rep(1,4),partition=x)})
l5 <- apply(sp,2, function(x){lp(evidence=a[5,],alpha=rep(1,4),partition=x)})

likelihood <- l1*l2*l3*l4*l5
likelihood <- likelihood/max(likelihood)
support <- log(likelihood)


###################################################
### code chunk number 14: plotlikelihoods
###################################################
par(mfcol=c(1,2))
plot(likelihood,ylab="likelihood",xlab="hypothesis")
abline(h = exp(-2))
plot(support,ylab="support",xlab="hypothesis")
abline(h = -2)


###################################################
### code chunk number 15: show.maximum.likelihood
###################################################
sp <- setparts(7)


###################################################
### code chunk number 16: setpartitions.Rnw:673-674
###################################################
sp[,which.max(support)]


###################################################
### code chunk number 17: likelihood.gt.2
###################################################
index <- order(support,decreasing=TRUE)
sp <- sp[,index]
support <- support[index]
dimnames(sp) <- list(
                     crime = paste("C",1:7, sep=""),
                     partition =  paste("H",1:ncol(sp),sep="")
                     )


###################################################
### code chunk number 18: setpartitions.Rnw:695-696
###################################################
 sp[, support > -2]


###################################################
### code chunk number 19: show.likelihoods
###################################################
support[support > -2]


