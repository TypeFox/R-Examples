rm(list = ls())

library(PPtree)
library(EMCluster)
library(MASS)
# source("./R/information.r")
# load("./data/crabs.em.rda")
load("./data/crabs.RndEM.rda")

x <- as.matrix(crabs[, 4:8])
# x <- t((t(x) - colMeans(x)) / sqrt(diag(var(x))))

da.s.all <- list()
for(k0 in 2:6){
  emobj <- ret.save[[k0]]
  var <- LTSigma2variance(emobj$LTSigma)
  x.pp <- x
  Sigma.pp <- var
  for(k.var in 1:dim(var)[3]){
    tmp <- eigen(var[,, k.var])
    Sigma.k.inv <- tmp$vector %*% diag(sqrt(1/tmp$values)) %*% t(tmp$vector)

    tmp <- x[emobj$class == k.var,]
    tmp.mu <- colMeans(tmp)
    tmp <- t(t(tmp) - tmp.mu) %*% Sigma.k.inv
    tmp <- t(t(tmp) + as.vector(tmp.mu))
    x.pp[emobj$class == k.var,] <- tmp

    # Sigma.pp[,, k.var] <- Sigma.k.inv
  }
  set.seed(1234)
  tmp.pp <- PP.optimize.random("LDA", 2, x.pp, emobj$class, std = FALSE)
  # plot(y %*% tmp.pp$proj.best, xlab = "PP1", ylab = "PP2",
  #      col = ret.save[[k0]]$class, pch = ret.save[[k0]]$class)
  x.new <- x %*% tmp.pp$proj.best
  mu.new <- emobj$Mu %*% tmp.pp$proj.best
  var.new <- array(0, dim = c(2, 2, dim(var)[3]))
  for(k.var in 1:dim(var)[3]){
    tmp <- t(tmp.pp$proj.best) %*% Sigma.pp[,, k.var] %*% tmp.pp$proj.best
    var.new[,, k.var] <- tmp
  }
  da.s <- list(pi = emobj$pi,
               Mu = mu.new,
               LTSigma = variance2LTSigma(var.new),
               class = emobj$class,
               nclass = emobj$nclass)
  da.s.all[[k0]] <- list(model = da.s, x = x.new)

  # color.class <- 1:11
  # postscript(file = paste("./plot/crabs_pp_k=", k0, ".eps", sep = ""),
  #            height = 6, width = 6, horizontal = FALSE)
  # plotem(da.s, x.new, xlab = "PP1", ylab = "PP2",
  #        main = paste("K = ", k0, sep = ""))
  # dev.off()
}

save(da.s.all, x, file = "data/crabs.pp.rda")
