plotNice2 <- function(cc2, d1, d2, XY="X", cex=0.7){
    
  plot(1, type="n", xlab=paste("loading ", d1, sep=""), ylab=paste("loading ", d2, sep=""), xlim=c(-1,1), ylim=c(-1,1))
  draw.circle(0, 0, 1, col="black")
  draw.circle(0, 0, 0.5, col="black")
  abline(v=0)
  abline(h=0)
  
  gidx <- c()
  if(XY=="X"){
    for(i in 1:nrow(cc2$corr.X.xscores)){
      x <- cc2$corr.X.xscores[i,d1]
      y <- cc2$corr.X.xscores[i,d2]
      
      rs <- x**2 + y**2
      
      if(rs >= 0.5**2){
        gidx <- c(gidx, i)        
        text(x=x, y=y, labels=as.character(rownames(cc2$corr.X.xscores)[i]), cex=cex, font=2)
      }
      if(rs < 0.5**2){
        points(x,y,pch=20, cex=cex, col="gray")
      }
    }  
  }
  if(XY=="Y"){
    for(i in 1:nrow(cc2$corr.Y.xscores)){
      x <- cc2$corr.Y.xscores[i,d1]
      y <- cc2$corr.Y.xscores[i,d2]
      
      rs <- x**2 + y**2
      
      if(rs >= 0.5**2){
        gidx <- c(gidx, i)
        text(x=x, y=y, labels=as.character(rownames(cc2$corr.Y.xscores)[i]), cex=cex, font=2)
      }
      if(rs < 0.5**2){
        points(x,y,pch=20, cex=cex, col="gray")
      }
    }  
  }
}
# ---
draw.circle <- function (x, y, r, col) {
    lines(cos(seq(0, 2 * pi, pi/180)) * r + x, sin(seq(0, 2 * 
        pi, pi/180)) * r + y, col = col)
}
# ---
plotNice2.indiv <- function(cc2, d1, d2,cex=0.7,cols=NULL, pchs=NULL){
  
  temp = range(cc2$xscores[,d1:d2])
  lims = c(floor(temp[1]), ceiling(temp[2]))
  
  plot(1, type="n", xlab=paste("Dimension ", d1, sep=""), ylab=paste("Dimension ", d2, sep=""), xlim=lims, ylim=lims)
  abline(v=0, lty=3)
  abline(h=0, lty=3)
  if(is.null(cols)){
    if(is.null(pchs)){
      text(cc2$xscores[, d1], cc2$xscores[, d2], rownames(cc2$xscores), cex=cex)}
    else{
      points(cc2$xscores[, d1], cc2$xscores[, d2], cex=cex, pch=pchs)}
  }
  else{
    if(is.null(pchs)){
      text(cc2$xscores[, d1], cc2$xscores[, d2], rownames(cc2$xscores), cex=cex, col=cols)}
    else{points(cc2$xscores[, d1], cc2$xscores[, d2], cex=cex, col=cols, pch=pchs)}
  }
}
# ----
mypamr.plotsurvival <- function (group, survival.time, censoring.status, cols, lty, lwd) 
{
  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
  junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
  pv <- 1 - pchisq(2 * (junk2$loglik[2] - junk2$loglik[1]), 
                   df = n.class - 1)
  plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd)
  legend(0.01 * max(survival.time), 0.2, col = cols, lty = lty, lwd=lwd, legend = as.character(1:n.class), box.lwd=1)
  text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
  return()
}
#