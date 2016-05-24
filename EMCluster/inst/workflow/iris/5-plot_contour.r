rm(list = ls())

library(EMCluster)
# source("./R/libContourPlotNew.r")
load("./data/iris.pp.rda")

for(k0 in 2:6){
  Pi <- da.s.all[[k0]]$model$pi
  Mu <- da.s.all[[k0]]$model$Mu
  S <- LTSigma2variance(da.s.all[[k0]]$model$LTSigma)
  class <- da.s.all[[k0]]$model$class
  da <- da.s.all[[k0]]$x

  pdf(paste("./plot/iris_ppcont_k=", k0, ".pdf", sep = ""),
      width = 6, height = 6)
  plotppcontour(da, Pi, Mu, S, class)
  dev.off()
}
