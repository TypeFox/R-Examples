### R code from vignette source 'Oncotree_vignette.Rnw'

###################################################
### code chunk number 1: Intro
###################################################
 options(width=100)
 ps.options(colormodel="rgb")


###################################################
### code chunk number 2: DataLoad
###################################################
 library(Oncotree)
 data(ov.cgh)
 str(ov.cgh)


###################################################
### code chunk number 3: TreeFit
###################################################
  ov.tree <- oncotree.fit(ov.cgh)


###################################################
### code chunk number 4: TreePrint
###################################################
  ov.tree
  plot(ov.tree, edge.weights="est")   


###################################################
### code chunk number 5: TreePrint2
###################################################
  pstree.oncotree(ov.tree, edge.weights="est", shape="oval")


###################################################
### code chunk number 6: TreePlot
###################################################
  plot(ov.tree, edge.weights="est")


###################################################
### code chunk number 7: TreeMarg
###################################################
  print(obs <- colMeans(ov.tree$data))
  print(est <- marginal.distr(ov.tree, with.errors=TRUE))
  #plot is in Figure 2
  barplot(rbind(obs[-1],est[-1]), beside=T,  legend.text=c("Observed","Fitted"),
	        main="Marginal frequencies of occurrence")


###################################################
### code chunk number 8: TreeMargPlot
###################################################
  barplot(rbind(obs[-1],est[-1]), beside=T,  legend.text=c("Observed","Fitted"),
	        main="Marginal frequencies of occurrence")


###################################################
### code chunk number 9: TreeDist
###################################################
  dd <- distribution.oncotree(ov.tree, with.errors=TRUE)
	head(dd)


###################################################
### code chunk number 10: Marg2way
###################################################
  #estimated probabilities of 2 events
  print(est2way <- t(data.matrix(dd[2:8])) %*% diag(dd$Prob) %*% data.matrix(dd[2:8]))
  #observed probabilities of 2 events
  print(obs2way <- t(ov.tree$data[,-1]) %*% ov.tree$data[,-1]/nrow(ov.tree$data))
  oe.diff <- obs2way-est2way
  oe.diff[lower.tri(oe.diff)] <- NA  #clear half of symmetric matrix for plotting
  
  require(lattice)  #the plot is in Figure 3
  levelplot(oe.diff, xlab="", ylab="", scales=list(x=list(alternating=2), tck=0),
	          main="Observed - Expected probabilities of co-occurrence of events")


###################################################
### code chunk number 11: Marg2wayPlot
###################################################
  print(levelplot(oe.diff, xlab="", ylab="", scales=list(x=list(alternating=2), tck=0),
	          main="Observed - Expected probabilities of co-occurrence of events"))


###################################################
### code chunk number 12: BootNP
###################################################
  set.seed(43636)
  ov.boot <- bootstrap.oncotree(ov.tree, type="nonparam", R=1000)
  ov.boot
  opar <- par(mfrow=c(3,2))   #the plot is in Figure 4
    plot(ov.boot, minfreq=45)
  par(opar)  


###################################################
### code chunk number 13: BootNPplot
###################################################
  opar <- par(mfrow=c(3,2))   #the plot is in Figure 4
    plot(ov.boot, minfreq=45, cex=1)
  par(opar)  


###################################################
### code chunk number 14: BootFreq
###################################################
  ov.boot$parent.freq


