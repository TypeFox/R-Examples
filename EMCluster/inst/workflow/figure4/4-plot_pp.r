rm(list = ls())

library(PPtree)
library(EMCluster)
# source("./R/information.r")
load("./data/simu.data.rda")
# load("./data/simu.em.rda")
load("./data/simu.RndEM.rda")

x <- da$X

da.s.all <- list()
for(k0 in 2:7){
  emobj <- ret.save[[k0]]
  da.s <- list(pi = emobj$pi,
               Mu = emobj$Mu,
               LTSigma = emobj$LTSigma,
               class = emobj$class,
               nclass = emobj$nclass)
  da.s.all[[k0]] <- list(model = da.s, x = x)

  color.class <- 1:11
  postscript(file = paste("./plot/figure4_pp_k=", k0, ".eps", sep = ""),
             height = 6, width = 6, horizontal = FALSE)
  plotem(da.s, x, xlab = "PP1", ylab = "PP2",
         main = paste("K = ", k0, sep = ""))
  dev.off()
}

save(da.s.all, x, file = "data/simu.pp.rda")
