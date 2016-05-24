#
# This demo illustrate a focusing behaviour of imprecise probabilities
#
# updated on 2015.11.27

set.seed(16979238)
lambda <- 1
n <- 2e1
y <- rpois(n=n, lambda=lambda)
y <- y[y>0]

pmat <- rbind(c(1.0,2.0),c(1.2,1.8),c(1.4,1.6),c(1.6,1.4),c(1.8,1.2),c(2.0,1.0))
rownames(pmat) <- paste("p", 1:nrow(pmat), sep="")
colnames(pmat) <- paste("d", 1:ncol(pmat), sep="")
op <- iprior(pmat=pmat)
op$vtx <- pmat
# plot(op)
op <- update(op, y=y[1], ztrunc=FALSE) # ztrunc=TRUE
# plot(op, after.obs=TRUE)
pbox(x=op, min=-5, max=3, main="ztrunc=FALSE")
sop <- summary(op)
sop

op <- update(op, y=y[1], ztrunc=TRUE)
sop <- summary(op)
sop
pbox(x=op, min=-5, max=3, main="ztrunc=TRUE")
