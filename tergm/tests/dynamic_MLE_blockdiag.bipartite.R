#  File tests/dynamic_MLE_blockdiag.bipartite.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(tergm)

tolerance<-0.05
n<-20
m<-13

theta<--.4

prop.weights <- c("default", "random")

logit<-function(p) log(p/(1-p))

block.dyadcount<-function(y, na.omit=TRUE){
  a <- y %v% "a"
  M <- outer(a,a,"==")
  M[1:m,1:m]<-0
  M[(m+1):n,(m+1):n]<-0
  M[lower.tri(M, TRUE)]<-0
  if(na.omit) M[as.edgelist(is.na(y))] <- 0
  sum(M)
}

form.mle<-function(y0,y1){
  logit(network.edgecount(y1-y0,na.omit=TRUE)/(block.dyadcount(y1)-network.edgecount(y0-is.na(y1))))
}

diss.mle<-function(y0,y1){
  -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
}

y0 <- network.initialize(n, directed=FALSE, bipartite=m)
a <- rep(1:6,1:6)[1:n]
a <- unlist(split(a, rep(1:2, n/2)))
a <- c(sort(a[1:m]), sort(a[-(1:m)]))
y0 %v% "a" <- a

set.seed(321)
y0<-standardize.network(simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2)))

cat("Complete data:\n")

set.seed(123)
y1<-standardize.network(simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2)))

# Force CMPLE
set.seed(543)
fit<-stergm(list(y0,y1), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMPLE", times=c(1,2))

stopifnot(fit$estimate=="CMPLE", fit$formation.fit$estimate=="MPLE", fit$dissolution.fit$estimate=="MPLE")
stopifnot(all.equal(form.mle(y0,y1), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Autodetected CMPLE
set.seed(543)
fit<-stergm(list(y0,y1), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE", times=c(1,2))

stopifnot(fit$estimate=="CMLE", is.null(fit$formation.fit$sample), is.null(fit$dissolution.fit$sample))
stopifnot(all.equal(form.mle(y0,y1), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")
set.seed(543)
fit<-stergm(list(y0,y1), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE", control=control.stergm(CMLE.control=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight)), times=c(1,2))

stopifnot(fit$estimate=="CMLE", fit$formation.fit$estimate=="MLE", fit$dissolution.fit$estimate=="MLE")
stopifnot(all.equal(form.mle(y0,y1), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))
}

cat("Missing data:\n")

y1m<-network.copy(y1)
set.seed(765)
e <- as.edgelist(y1)[1,]
y1m[e[1], e[2]] <- NA
y1m[1,m+2] <- NA

# Force CMPLE
set.seed(765)
fit<-stergm(list(y0,y1m), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMPLE", times=c(1,2))

stopifnot(fit$estimate=="CMPLE", fit$formation.fit$estimate=="MPLE", fit$dissolution.fit$estimate=="MPLE")
stopifnot(all.equal(form.mle(y0,y1m), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1m), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Autodetected CMPLE
set.seed(765)
fit<-stergm(list(y0,y1m), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE", times=c(1,2))

stopifnot(fit$estimate=="CMLE", is.null(fit$formation.fit$sample), is.null(fit$dissolution.fit$sample))
stopifnot(all.equal(form.mle(y0,y1m), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1m), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")
set.seed(234)
fit<-stergm(list(y0,y1m), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE", control=control.stergm(CMLE.control=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight)), times=c(1,2))

stopifnot(fit$estimate=="CMLE", fit$formation.fit$estimate=="MLE", fit$dissolution.fit$estimate=="MLE")
stopifnot(all.equal(form.mle(y0,y1m), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1m), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))
}

}, "dynamic MLE with block-diagonal constraints")
