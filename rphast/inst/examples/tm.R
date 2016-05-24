tree <- "((human:0.01, chimp:0.01):0.03, mouse:0.3)"
subst.mod <- "JC69"
rate.mat <- matrix(runif(16), nrow=4, ncol=4)
for (i in 1:4)
  rate.mat[i,i] <- -sum(rate.mat[i,-i])
backgd <- runif(4)
backgd <- backgd/sum(backgd)
alphabet <- "ACGT"
t <- tm(tree, subst.mod, rate.mat, backgd, alphabet)
t

nratecats <- 3
alpha <- 1.5
rate.consts <- runif(nratecats, max=3.0)
root.leaf <- "human"
t <- tm(tree, subst.mod, rate.matrix=rate.mat,
        backgd=backgd, alphabet=alphabet,
        nratecats=nratecats, alpha=alpha,
        rate.consts=rate.consts, root.leaf=root.leaf)
t

