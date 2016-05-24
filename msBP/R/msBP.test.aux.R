msBP.posteriorH0 <- function(i, group0, group1, pool, priorH0, a, b, plot.res=FALSE)
{
priorH1 <- 1 - priorH0
odd <- priorH1/priorH0
n0 <- group0$mcmcsamples$n[i,]
n1 <- group1$mcmcsamples$n[i,]
n <-  pool$mcmcsamples$n[i,]
r0 <- group0$mcmcsamples$r[i,]
r1 <- group1$mcmcsamples$r[i,]
r <-  pool$mcmcsamples$r[i,]
v0 <- group0$mcmcsamples$v[i,]
v1 <- group1$mcmcsamples$v[i,]
v <-  pool$mcmcsamples$v[i,]
maxScale <- group0$postmean$n$max.s 

msBP.nesting(n,r,v,n0,r0,v0,n1,r1,v1,priorH0,a,b,maxScale)
}
#-------------------------------------------------------------------------------
msBP.clus2prob <- function(group0, group1, pool, priorH0, a, b, plot.res=FALSE)
{
priorH1 <- 1 - priorH0
odd <- priorH1/priorH0

n0 <- tree2vec(group0$n)
n1 <- tree2vec(group1$n)
n <-  tree2vec(pool$n)
r0 <- tree2vec(group0$r)
r1 <- tree2vec(group1$r)
r <-  tree2vec(pool$r)
v0 <- tree2vec(group0$v)
v1 <- tree2vec(group1$v)
v <-  tree2vec(pool$v)
maxScale <- group0$n$max.s 

l0 <- v0-n0-r0
l1 <- v1-n1-r1
l <- v-n-r

msBP.nesting(n,r,v,n0,r0,v0,n1,r1,v1,priorH0,a,b,maxScale)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
msBP.nesting <- function(n,r,v,n0,r0,v0,n1,r1,v1,priorH0,a,b,maxScale)
{
priorH1 <- 1 - priorH0
odd <- priorH1/priorH0
l0 <- v0-n0-r0
l1 <- v1-n1-r1
l <- v-n-r

factor0_A <- lgamma(n0 + 1)+lgamma(a+v0-n0)-lgamma(a+v0+1)
factor1_A <- lgamma(n1 + 1)+lgamma(a+v1-n1)-lgamma(a+v1+1)
factor_A <-  lgamma(n  + 1)+lgamma(a+v -n )-lgamma(a+v +1)

factor0_B <- lgamma(b+r0)+lgamma(b+l0)-lgamma(2*b+v0-n0)
factor1_B <- lgamma(b+r1)+lgamma(b+l1)-lgamma(2*b+v1-n1)
factor_B <-  lgamma(b+r )+lgamma(b+l )-lgamma(2*b+v -n )

nestH1_0 <- sapply(vec2tree(factor0_A)$T,sum, na.rm=TRUE) + sapply(vec2tree(factor0_B)$T,sum, na.rm=TRUE)
nestH1_1 <- sapply(vec2tree(factor1_A)$T,sum, na.rm=TRUE) + sapply(vec2tree(factor1_B)$T,sum, na.rm=TRUE)
nestH000 <- sapply(vec2tree(factor_A)$T,sum, na.rm=TRUE)  + sapply(vec2tree(factor_B)$T,sum, na.rm=TRUE)

nestH1_0 <- nestH1_0[1:maxScale]
nestH1_1 <- nestH1_1[1:maxScale]
nestH000 <- nestH000[1:maxScale]

nestH1 <-    ((2^(1:maxScale-1)) * (lgamma(a+1) - lgamma(a) + lgamma(2*b) - 2*lgamma(b))) + nestH1_0 + nestH1_1 
nestH0 <-    ((2^(1:maxScale-1)) * (lgamma(a+1) - lgamma(a) + lgamma(2*b) - 2*lgamma(b))) + nestH000

postH0 <- 1/(1 + odd*exp(nestH1-nestH0))
BF <- exp(nestH0 -  nestH1)
rbind(postH0=postH0[1:(maxScale)], nestH0=nestH0[1:(maxScale)], nestH1 = nestH1[1:(maxScale)], BF=BF[1:(maxScale)])
}
#---
