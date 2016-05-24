### R code from vignette source 'surface_tutorial.Rnw'

###################################################
### code chunk number 1: surface_tutorial.Rnw:23-28
###################################################
library(surface)

data(surfaceDemo)
tree<-surfaceDemo$tree
dat<-surfaceDemo$sim$dat


###################################################
### code chunk number 2: surface_tutorial.Rnw:32-35
###################################################
tree<-nameNodes(tree)
olist<-convertTreeData(tree,dat)
otree<-olist[[1]]; odata<-olist[[2]]


###################################################
### code chunk number 3: surface_tutorial.Rnw:43-46
###################################################
fwd<-surfaceForward(otree, odata, aic_threshold = 0, exclude = 0, 
verbose = FALSE, plotaic = FALSE)
k<-length(fwd)


###################################################
### code chunk number 4: surface_tutorial.Rnw:49-52
###################################################
fsum<-surfaceSummary(fwd)
names(fsum)
fsum$aics


###################################################
### code chunk number 5: surface_tutorial.Rnw:56-60
###################################################
bwd<-surfaceBackward(otree, odata, starting_model = fwd[[k]], aic_threshold = 0, 
only_best = TRUE, verbose = FALSE, plotaic = FALSE)
bsum<-surfaceSummary(bwd)
kk<-length(bwd)


###################################################
### code chunk number 6: surface_tutorial.Rnw:63-67
###################################################
bsum$alpha
bsum$sigma_squared
bsum$theta
bsum$n_regimes


###################################################
### code chunk number 7: surface_tutorial.Rnw:74-75
###################################################
surfaceTreePlot(tree, bwd[[kk]], labelshifts = T)


###################################################
### code chunk number 8: surface_tutorial.Rnw:77-80
###################################################
par(mfrow=c(1,2), mai=c(0.8,0.8,0.2,0.2))
surfaceTraitPlot(dat, bwd[[kk]], whattraits = c(1,2))
surfaceTraitPlot(dat, bwd[[kk]], whattraits = c(3,2))


###################################################
### code chunk number 9: surface_tutorial.Rnw:86-88
###################################################
bm<-startingModel(otree,odata,brownian=TRUE)
ou1<-startingModel(otree,odata)


###################################################
### code chunk number 10: surface_tutorial.Rnw:91-92
###################################################
H12<-startingModel(otree,odata,shifts=c("26"="H1","13"="H1","5"="H2","19"="H2"))


###################################################
### code chunk number 11: surface_tutorial.Rnw:97-101
###################################################
surfaceAICPlot(fwd, bwd)
abline(h=bm[[1]]$aic,lty="longdash")
abline(h=H12[[1]]$aic,lty="longdash")
text(c(6,6),c(bm[[1]]$aic, ou1[[1]]$aic, H12[[1]]$aic)-2,c("BM","OU1","H12"),cex=0.5) 


###################################################
### code chunk number 12: surface_tutorial.Rnw:107-110
###################################################
truefit<-surfaceDemo$sim$fit
propRegMatch(truefit, bwd[[kk]]$fit, internal = FALSE)  
propRegMatch(truefit, bwd[[kk]]$fit, internal = TRUE)  


###################################################
### code chunk number 13: surface_tutorial.Rnw:114-117
###################################################
set.seed(10)
newsim<-surfaceSimulate(tree, type="hansen-fit", hansenfit=fwd[[k]]$fit, 
shifts=fwd[[k]]$savedshifts, sample_optima=TRUE)


###################################################
### code chunk number 14: surface_tutorial.Rnw:121-124
###################################################
par(mfrow=c(1,2),mai=c(0.8,0.8,0.2,0.2))
surfaceTraitPlot(newsim$data, newsim, whattraits = c(1,2), convcol = FALSE)
surfaceTraitPlot(newsim$data, newsim, whattraits = c(3,2), convcol = FALSE)


###################################################
### code chunk number 15: surface_tutorial.Rnw:131-136
###################################################
newout<-runSurface(tree, newsim$dat, only_best = TRUE)
newsum<-surfaceSummary(newout$bwd)
newkk<-length(newout$bwd)
newsum$n_regimes
bsum$n_regimes


###################################################
### code chunk number 16: surface_tutorial.Rnw:139-142
###################################################
par(mfrow=c(1,2),mai=c(0.8,0.8,0.2,0.2))
surfaceTraitPlot(newsim$data, newout$bwd[[newkk]], whattraits = c(1,2))
surfaceTraitPlot(newsim$data, newout$bwd[[newkk]], whattraits = c(3,2))


