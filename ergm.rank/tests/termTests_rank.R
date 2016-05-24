#  File tests/termTests_rank.R in package ergm.rank, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2016 Statnet Commons
#######################################################################
library(ergm.rank)

local1.nonconformity <- function(m){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(k) sapply(1:n, function(l) if(length(unique(c(i,j,k,l)))==4) (m[i,l]>m[i,j])*(m[l,j]>m[l,k])*(1-(m[i,j]>m[i,k])) else 0 )))))
}
local2.nonconformity <- function(m){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(k) sapply(1:n, function(l) if(length(unique(c(i,j,k,l)))==4) (m[i,l]>m[i,k])*(m[l,j]>m[l,k])*(1-(m[i,j]>m[i,k])) else 0 )))))
}
localAND.nonconformity <- function(m){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(k) sapply(1:n, function(l) if(length(unique(c(i,j,k,l)))==4) (m[i,l]>m[i,k] && m[i,l]>m[i,j])*(m[l,j]>m[l,k])*(1-(m[i,j]>m[i,k])) else 0 )))))
}
global.nonconformity <- function(m){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(k) sapply(1:n, function(l) if(length(unique(c(i,j,k,l)))==4) (m[l,j]>m[l,k])*(1-(m[i,j]>m[i,k])) else 0 )))))
}
deference <- function(m){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(l) if(length(unique(c(i,j,l)))==3) (m[l,j]>m[l,i])*(m[i,l]>m[i,j]) else 0 ))))
}
nodeicov <- function(m, x){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(k) if(length(unique(c(i,j,k)))==3) (m[i,j]>m[i,k])*(x[j]-x[k]) else 0 ))))
}
edgecov <- function(m, x){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(k) if(length(unique(c(i,j,k)))==3) (m[i,j]>m[i,k])*(x[i,j]-x[i,k]) else 0 ))))
}
inconsistency <- function(m, x, w=array(1,c(nrow(m),ncol(m),ncol(m)))){
  n <- nrow(m)
  sum(sapply(1:n, function(i) sapply(1:n, function(j) sapply(1:n, function(k) if(length(unique(c(i,j,k)))==3) ((m[i,j]>m[i,k])*(x[i,j]<=x[i,k])+(m[i,j]<=m[i,k])*(x[i,j]>x[i,k]))*(w[i,j,k]) else 0 ))))
}
n <- 7
S <- 10

nw0 <- network.initialize(n,dir=TRUE)
nw0[,, names.eval="r", add.edges=TRUE] <- m0 <- matrix(sample.int(n, n*n, replace=TRUE), n, n)
nw0 %v% "v" <- xv <- rnorm(n)
nw0 %n% "m" <- xm <- matrix(rnorm(n*n),n,n)
nw0 %n% "a" <- xa <- array(rnorm(n*n*n),c(n,n,n))

nws <- simulate(nw0~rank.nonconformity("all")+
                  rank.nonconformity("local1")+
                  rank.nonconformity("local2")+
                  rank.nonconformity("localAND")+
                  rank.deference+
                  rank.nodeicov("v")+
                  rank.edgecov("m")+
                  rank.inconsistency(nw0,"r",xa),
                coef=rep(0,8),response="r", reference=~DiscUnif(1, n-1), nsim=S, statsonly=FALSE)

d.stats <- attr(nws,"stats")

c.stats <- t(sapply(lapply(nws,as.matrix,attrname="r"),
                    function(m) c(global.nonconformity(m),
                                  local1.nonconformity(m),
                                  local2.nonconformity(m),
                                  localAND.nonconformity(m),
                                  deference(m),
                                  nodeicov(m, xv),
                                  edgecov(m, xm),
                                  inconsistency(m, m0, xa))))

s.stats <- summary(nws~rank.nonconformity("all")+
                  rank.nonconformity("local1")+
                  rank.nonconformity("local2")+
                  rank.nonconformity("localAND")+
                  rank.deference+
                  rank.nodeicov("v")+
                  rank.edgecov("m")+
                  rank.inconsistency(nw0,"r",xa), response="r")

stopifnot(all.equal(d.stats,c.stats,check.attributes=FALSE),
          all.equal(d.stats,s.stats,check.attributes=FALSE),
          all.equal(s.stats,c.stats,check.attributes=FALSE))

nw1 <- nws[[length(nws)]]

nws <- simulate(nw1~rank.nonconformity("all")+
                  rank.nonconformity("local1")+
                  rank.nonconformity("local2")+
                  rank.nonconformity("localAND")+
                  rank.deference+
                  rank.nodeicov("v")+
                  rank.edgecov("m")+
                  rank.inconsistency(nw0,"r",xa), coef=rep(0,8),response="r", reference=~CompleteOrder, nsim=S, statsonly=FALSE)

d.stats <- attr(nws,"stats")

c.stats <- t(sapply(lapply(nws,as.matrix,attrname="r"),
                    function(m) c(global.nonconformity(m),
                                  local1.nonconformity(m),
                                  local2.nonconformity(m),
                                  localAND.nonconformity(m),
                                  deference(m),
                                  nodeicov(m, xv),
                                  edgecov(m, xm),
                                  inconsistency(m, m0, xa))))

s.stats <- summary(nws~rank.nonconformity("all")+
                  rank.nonconformity("local1")+
                  rank.nonconformity("local2")+
                  rank.nonconformity("localAND")+
                  rank.deference+
                  rank.nodeicov("v")+
                  rank.edgecov("m")+
                  rank.inconsistency(nw0,"r",xa), response="r")

stopifnot(all.equal(d.stats,c.stats,check.attributes=FALSE),
          all.equal(d.stats,s.stats,check.attributes=FALSE),
          all.equal(s.stats,c.stats,check.attributes=FALSE))
