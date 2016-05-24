
source("lungLoad.r")
source("cvKfold.r")

########################################################################################################################
## Boston cohort

npars <- 1
ncats <- 3
pathwayFile <- "c2.cp.kegg.v3.0.symbols.gmt"

avrres <- NULL
cvmean <- NULL
cvstderr <- NULL
cvpathway <- NULL
for(nn in (1+length(avrres)):186) {
  d <- load.gsea(pathwayFile,"Boston",nn,npars,ncats)
  if(nrow(d$cdata)*ncol(d$cdata) <= 0 || nrow(d$cdata) <= 1)
    next
  ddata <- cnDiscretize(d$cdata, ncats, "uniform")
  res10 <- NULL
  for(k in 1:10)
    res10 <- rbind(res10, classification.cv10fold(d$cdata,d$cls,d$pathway,ddata,npars,d$nodePars,ncats,d$nodeCats,d$nodeOrder,ker="radial",0.03,folds=10,echo=FALSE))
  sc <- max(mean(unlist(res10[,1])), mean(unlist(res10[,2])), mean(unlist(res10[,3])), mean(unlist(res10[,4])), mean(unlist(res10[,5])), mean(unlist(res10[,6])), mean(unlist(res10[,7])))
  if(sc > 60) {
    for(k in 1:10)
      res10 <- rbind(res10, classification.cv10fold(d$cdata,d$cls,d$pathway,ddata,npars,d$nodePars,ncats,d$nodeCats,d$nodeOrder,ker="radial",0.03,folds=10,echo=FALSE))
    err <- try(postscript(paste("./kegg/boston_", d$pathname, "_p1c3.ps", sep="")), silent=TRUE)
    if(is.null(err)) {
      boxplot(list(NB=unlist(res10[,1]), BN.0=unlist(res10[,2]), BN.BIC=unlist(res10[,3]), BN.AIC=unlist(res10[,4]), SVM=unlist(res10[,5]), LASSO=unlist(res10[,6]), LAR=unlist(res10[,7])), main=paste("Lung, Boston, ", d$pathname, sep=""))
      dev.off()
    }
  }
  cvmean <- rbind(cvmean, c(mean(unlist(res10[,1])), mean(unlist(res10[,2])), mean(unlist(res10[,3])), mean(unlist(res10[,4])), mean(unlist(res10[,5])), mean(unlist(res10[,6])), mean(unlist(res10[,7]))))
  cvstderr <- rbind(cvstderr, c(sqrt(var(unlist(res10[,1]))), sqrt(var(unlist(res10[,2]))), sqrt(var(unlist(res10[,3]))), sqrt(var(unlist(res10[,4]))), sqrt(var(unlist(res10[,5]))), sqrt(var(unlist(res10[,6]))), sqrt(var(unlist(res10[,7])))))
  sc <- max(mean(unlist(res10[,1])), mean(unlist(res10[,2])), mean(unlist(res10[,3])), mean(unlist(res10[,4])), mean(unlist(res10[,5])), mean(unlist(res10[,6])), mean(unlist(res10[,7])))
  avrres <- c(avrres, sc)
  cvpathway <- c(cvpathway, d$pathname)
  cat(nn,": ", sc,"\n")
}

tt <- cbind(cvpathway,cvmean,cvstderr, avrres)
colnames(tt) <- c("pathway", "NB", "BN.0", "BN.BIC", "BN.AIC", "SVM", "LASSO", "LAR", "STD.NB", "STD.BN.0", "STD.BN.BIC", "STD.BN.AIC", "STD.SVM", "STD.LASSO", "STD.LAR", "maxscore")
write.csv(tt, "boston_cv10_p1c3.csv")

########################################################################################################################
## Boston and Michigan combined 

avrres <- NULL
cvmean <- NULL
cvstderr <- NULL
cvpathway <- NULL
for(nn in (1+length(avrres)):186) {
  bd <- load.gsea(pathwayFile,"Boston",nn,npars,ncats)
  md <- load.gsea(pathwayFile,"Michigan",nn, npars, ncats)
  md$cls <- md$cls[1:48]
  md$cdata <- md$cdata[,1:48]
  pathway <- bd$pathway
  rgen <- NULL
  for(gen in 1:length(pathway))
    if(length(which(md$pathway == pathway[gen]))<1)
      rgen <- c(rgen, gen)
  if(!is.null(rgen)) {
    pathway <- pathway[-rgen]
    bd$cdata <- bd$cdata[pathway,]
    md$cdata <- md$cdata[pathway,]
    bd$nodeOrder <- bd$nodeOrder[pathway]
    md$nodeOrder <- md$nodeOrder[pathway]
    bd$pathway <- pathway
    md$pathway <- pathway
  }
  names(bd$nodeOrder) <- pathway
  names(md$nodeOrder) <- pathway
  clslevs <- range(bd$cls)
  numnodes <- length(pathway)
  bd$cdata <- bd$cdata - apply(bd$cdata, 1, "mean")
  bd$cdata <- bd$cdata/sqrt(apply(bd$cdata, 1, "var")/ncol(bd$cdata))
  md$cdata <- md$cdata - apply(md$cdata, 1, "mean")
  md$cdata <- md$cdata/sqrt(apply(md$cdata, 1, "var")/ncol(md$cdata))
  ## combine
  d <- bd
  d$cdata <- cbind(bd$cdata, md$cdata)
  d$cls <- c(bd$cls,md$cls)
  if(nrow(d$cdata)*ncol(d$cdata) <= 0 || nrow(d$cdata) <= 1)
    next
  ddata <- cnDiscretize(d$cdata, ncats, "uniform")
  res10 <- NULL
  for(k in 1:5)
    res10 <- rbind(res10, classification.cv10fold(d$cdata,d$cls,d$pathway,ddata,npars,d$nodePars,ncats,d$nodeCats,d$nodeOrder, ker="radial",0.03, folds=10,echo=TRUE))
  sc <- max(mean(unlist(res10[,1])), mean(unlist(res10[,2])), mean(unlist(res10[,3])), mean(unlist(res10[,4])), mean(unlist(res10[,5])), mean(unlist(res10[,6])), mean(unlist(res10[,7])))
  if(sc > 60) {
    for(k in 1:10)
      res10 <- rbind(res10, classification.cv10fold(d$cdata,d$cls,d$pathway,ddata,npars,d$nodePars,ncats,d$nodeCats,d$nodeOrder, ker="radial",0.03, folds=10,echo=TRUE))
    err <- try(postscript(paste("./kegg/bm_", d$pathname, "_p1c3.ps", sep="")), silent=TRUE)
    if(is.null(err)) {
      boxplot(list(NB=unlist(res10[,1]), BN.0=unlist(res10[,2]), BN.BIC=unlist(res10[,3]), BN.AIC=unlist(res10[,4]), SVM=unlist(res10[,5]), LASSO=unlist(res10[,6]), LAR=unlist(res10[,7])), main=paste("Lung, ", d$pathname, sep=""))
      dev.off()
    }
  }
  cvmean <- rbind(cvmean, c(mean(unlist(res10[,1])), mean(unlist(res10[,2])), mean(unlist(res10[,3])), mean(unlist(res10[,4])), mean(unlist(res10[,5])), mean(unlist(res10[,6])), mean(unlist(res10[,7]))))
  cvstderr <- rbind(cvstderr, c(sqrt(var(unlist(res10[,1]))), sqrt(var(unlist(res10[,2]))), sqrt(var(unlist(res10[,3]))), sqrt(var(unlist(res10[,4]))), sqrt(var(unlist(res10[,5]))), sqrt(var(unlist(res10[,6]))), sqrt(var(unlist(res10[,7])))))
  sc <- max(mean(unlist(res10[,1])), mean(unlist(res10[,2])), mean(unlist(res10[,3])), mean(unlist(res10[,4])), mean(unlist(res10[,5])), mean(unlist(res10[,6])), mean(unlist(res10[,7])))
  avrres <- c(avrres, sc)
  cvpathway <- c(cvpathway, d$pathname)
  cat(nn,": ", sc,"\n")
}

tt <- cbind(cvpathway,cvmean,cvstderr, avrres)
colnames(tt) <- c("pathway", "NB", "NET.NULL", "NET.BIC", "NET.AIC", "SVM", "LASSO", "LAR", "STD.NB", "STD.NET.NULL", "STD.NET.BIC", "STD.NET.AIC", "STD.SVM", "STD.LASSO", "STD.LAR", "maxscore")
write.csv(tt, "bm_cv10_p1c3.csv")


