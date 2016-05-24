# 
# A mixture of two populations
# 
# updated on 2015.11.29

library(dsample)
expr <- expression((x1*(1-x2))^5 * (x2*(1-x1))^3 * (1-x1*(1-x2)-x2*(1-x1))^37)
sets <- list(x1=runif(1e5), x2=runif(1e5))
smp <- dsample(expr=expr, rpmat=sets, nk=1e4, n=1e3)
op <- summary(smp, n=10, k=2)
op$means
op$modes
do.call(cbind, lapply(split(op$X, op$grp), colMeans))
plot(op)

