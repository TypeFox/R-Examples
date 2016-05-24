##############################################################################
plot.cit <- function(sim.out, alpha, my.cex, main="")
{
  par(mar=c(5,4.5,4,4)+0.1)
  r2s <- sim.out$R2s
  plot(r2s[,1],r2s[,2],xlab=expression(paste(R^{2}*(list(Y[1],Q)))),
    ylab=expression(paste(R^{2}*(list(Y[2],Q)))),cex=my.cex,col="grey",
    lwd=2,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=main)
  abline(a=0,b=1,lwd=2)
  pvals <- sim.out$pval.cit
  #discoveries.G <- which((pvals[,1] < alpha) & (pvals[,2] < alpha))
  discoveries.M1 <- which((pvals[,1] < alpha) & (pvals[,2] > alpha))
  discoveries.M2 <- which((pvals[,1] > alpha) & (pvals[,2] < alpha))
  discoveries.M3 <- which((pvals[,1] > alpha) & (pvals[,2] > alpha))
  n <- nrow(r2s)
  aux <- rep(NA,n)
  aux[discoveries.M1] <- 1
  aux[discoveries.M2] <- 2
  aux[discoveries.M3] <- 3
  aux2 <- which(!is.na(aux))
  mycols <- c("blue","red","yellow")
  for(i in aux2)
    points(r2s[i,1],r2s[i,2],col=mycols[aux[i]],cex=my.cex,lwd=2)
}
##############################################################################
plot.bic.aic <- function(sim.out, my.cex, main="", penalty="BIC")
{
  par(mar=c(5,4.5,4,4)+0.1)
  r2s <- sim.out$R2s
  plot(r2s[,1],r2s[,2],xlab=expression(paste(R^{2}*(list(Y[1],Q)))),
    ylab=expression(paste(R^{2}*(list(Y[2],Q)))),cex=my.cex,type="n",
    lwd=2,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=main)
  abline(a=0,b=1,lwd=2)
  if (penalty == "BIC")
    ICs <- sim.out$BICs
  if (penalty == "AIC")
    ICs <- sim.out$AICs
  rank.ICs <- t(apply(ICs,1,rank))
  discoveries.M1 <- which(rank.ICs[,1]==1)
  discoveries.M2 <- which(rank.ICs[,2]==1)
  discoveries.M3 <- which(rank.ICs[,3]==1)
  discoveries.M4 <- which(rank.ICs[,4]==1)
  n <- nrow(r2s)
  aux <- rep(NA,n)
  aux[discoveries.M1] <- 1
  aux[discoveries.M2] <- 2
  aux[discoveries.M3] <- 3
  aux[discoveries.M4] <- 4
  mycols <- c("blue","red","green","black")
  for(i in 1:n)
    points(r2s[i,1],r2s[i,2],col=mycols[aux[i]],cex=my.cex,lwd=2)
}
##############################################################################
plot.par.cmst.joint <- function(sim.out, alpha, my.cex, 
  main="", penalty=c("BIC","AIC"))
{
  par(mar=c(5,4.5,4,4)+0.1)
  r2s <- sim.out$R2s
  plot(r2s[,1],r2s[,2],xlab=expression(paste(R^{2}*(list(Y[1],Q)))),
    ylab=expression(paste(R^{2}*(list(Y[2],Q)))),cex=my.cex,col="grey",
    lwd=2,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=main)
  abline(a=0,b=1,lwd=2)
  if(penalty=="BIC")
    pvals <- sim.out$pval.par.cmst.joint.BIC
  if(penalty=="AIC")
    pvals <- sim.out$pval.par.cmst.joint.AIC
  discoveries.M1 <- which(pvals[,1] <= alpha)
  discoveries.M2 <- which(pvals[,2] <= alpha)
  discoveries.M3 <- which(pvals[,3] <= alpha)
  discoveries.M4 <- which(pvals[,4] <= alpha)
  n <- nrow(r2s)
  aux <- rep(NA,n)
  aux[discoveries.M1] <- 1
  aux[discoveries.M2] <- 2
  aux[discoveries.M3] <- 3
  aux[discoveries.M4] <- 4
  mycols <- c("blue","red","green","black")
  aux2 <- which(!is.na(aux))
  for(i in aux2)
    points(r2s[i,1],r2s[i,2],col=mycols[aux[i]],cex=my.cex,lwd=2)
}
##############################################################################
plot.par.cmst <- function(sim.out, alpha, my.cex, 
  main="", penalty=c("BIC","AIC"))
{
  par(mar=c(5,4.5,4,4)+0.1)
  r2s <- sim.out$R2s
  plot(r2s[,1],r2s[,2],xlab=expression(paste(R^{2}*(list(Y[1],Q)))),
    ylab=expression(paste(R^{2}*(list(Y[2],Q)))),cex=my.cex,col="grey",
    lwd=2,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=main)
  abline(a=0,b=1,lwd=2)
  if(penalty=="BIC")
    pvals <- sim.out$pval.par.cmst.iu.BIC
  if(penalty=="AIC")
    pvals <- sim.out$pval.par.cmst.iu.AIC
  discoveries.M1 <- which(pvals[,1] <= alpha)
  discoveries.M2 <- which(pvals[,2] <= alpha)
  discoveries.M3 <- which(pvals[,3] <= alpha)
  discoveries.M4 <- which(pvals[,4] <= alpha)
  n <- nrow(r2s)
  aux <- rep(NA,n)
  aux[discoveries.M1] <- 1
  aux[discoveries.M2] <- 2
  aux[discoveries.M3] <- 3
  aux[discoveries.M4] <- 4
  mycols <- c("blue","red","green","black")
  aux2 <- which(!is.na(aux))
  for(i in aux2)
    points(r2s[i,1],r2s[i,2],col=mycols[aux[i]],cex=my.cex,lwd=2)
}
##############################################################################
plot.non.par.cmst <- function(sim.out, alpha, my.cex, 
  main="", penalty=c("BIC","AIC"))
{
  par(mar=c(5,4.5,4,4)+0.1)
  r2s <- sim.out$R2s
  plot(r2s[,1],r2s[,2],xlab=expression(paste(R^{2}*(list(Y[1],Q)))),
    ylab=expression(paste(R^{2}*(list(Y[2],Q)))),cex=my.cex,col="grey",
    lwd=2,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=main)
  abline(a=0,b=1,lwd=2)
  if(penalty=="BIC")
    pvals <- sim.out$pval.non.par.cmst.iu.BIC
  if(penalty=="AIC")
    pvals <- sim.out$pval.non.par.cmst.iu.AIC
  discoveries.M1 <- which(pvals[,1] <= alpha)
  discoveries.M2 <- which(pvals[,2] <= alpha)
  discoveries.M3 <- which(pvals[,3] <= alpha)
  discoveries.M4 <- which(pvals[,4] <= alpha)
  n <- nrow(r2s)
  aux <- rep(NA,n)
  aux[discoveries.M1] <- 1
  aux[discoveries.M2] <- 2
  aux[discoveries.M3] <- 3
  aux[discoveries.M4] <- 4
  mycols <- c("blue","red","green","black")
  aux2 <- which(!is.na(aux))
  for(i in aux2)
    points(r2s[i,1],r2s[i,2],col=mycols[aux[i]],cex=my.cex,lwd=2)
}
