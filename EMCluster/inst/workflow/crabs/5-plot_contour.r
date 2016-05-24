rm(list = ls())

library(EMCluster)
# source("./R/libContourPlotNew.r")
load("./data/crabs.pp.rda")

for(k0 in 2:6){
  Pi <- da.s.all[[k0]]$model$pi
  Mu <- da.s.all[[k0]]$model$Mu
  S <- LTSigma2variance(da.s.all[[k0]]$model$LTSigma)
  class <- da.s.all[[k0]]$model$class
  da <- da.s.all[[k0]]$x

  # crabs has weird projection and needs a rotation to get a better view.
  if(k0 == 2) angle <- 40 / 180 * pi
  if(k0 == 3) angle <- 10 / 180 * pi
  if(k0 == 4) angle <- 45 / 180 * pi
  if(k0 == 5) angle <- -20 / 180 * pi
  if(k0 == 6) angle <- 0 / 180 * pi

  pdf(paste("./plot/crabs_ppcont_k=", k0, ".pdf", sep = ""),
      width = 6, height = 6)
  plotppcontour(da, Pi, Mu, S, class, angle = angle)
  dev.off()
}
