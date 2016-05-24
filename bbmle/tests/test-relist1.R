library(bbmle)
set.seed(1001)
f <- factor(rep(1:3,each=50))
kvals <- c(1,2,5)
muvals <- c(10,2,5)
y <- rnbinom(length(f),size=kvals[f],mu=muvals[f])
plot(y)

NLL <- function(p) {
  kvals <- p[1:3]
  muvals <- p[4:6]
  -sum(dnbinom(y,size=kvals[f],mu=muvals[f],log=TRUE))
}
parnames(NLL) <- c("k1","k2","k3","mu1","mu2","mu3")
svec <- c(kvals,muvals)
names(svec) <- parnames(NLL)
m1 <- mle2(NLL,start=svec,vecpar=TRUE)
