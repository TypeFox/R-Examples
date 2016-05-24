
plotFRBvars <- function(x, cumul=2, confmethod = c("BCA","basic"), npcs = min(10, length(x$eigval)))
{
FRBres <- x
confmethod <- match.arg(confmethod)
conf <- FRBres$conf
if (!(cumul %in% c(0,1,2))) stop("argument 'cumul' should be 0, 1 or 2")

p <- length(FRBres$eigval)
if (npcs>p) stop("argument npcs should be at most the number of variables")

if (cumul==2) par(mfrow=c(2,1))
else par(mfrow=c(1,1))

if ((cumul==1)|(cumul==2)) {
  Fvars <- FRBres$pvar
  if (confmethod=="basic") {
    FvarsLow <- FRBres$pvar.CI.basic[,1]
    FvarsHigh <- FRBres$pvar.CI.basic[,2]
  } else {
    FvarsLow <- FRBres$pvar.CI.bca[,1]
    FvarsHigh <- FRBres$pvar.CI.bca[,2]
  }
  Fvars[p] <- 1
  FvarsLow[p] <- 1
  FvarsHigh[p] <- 1

  vars <- Fvars[1:npcs]*100
  varsLow <- FvarsLow[1:npcs]*100
  varsHigh <- FvarsHigh[1:npcs]*100

  xlimits <- c(0.5,npcs+0.5)
  emptyrange <- min(varsLow)
  if (emptyrange < 33)
    ylimits <- c(0,105)
  else
    ylimits <- c(100 - (100-emptyrange)*1.5, 100 + (100-emptyrange)*.05)
  
  if (cumul==2) {
    ylab=""
    main=paste("Cumulative % of variance (+ ",conf*100,"% ",confmethod," confidence limits)",sep="")
  }
  else {
    ylab="Cumulative % of variance"
    main=paste("Cumulative % of variance (+ ",conf*100,"% ",confmethod," confidence limits)",sep="")
  }
  plot(1:npcs, vars, pch=20, type="b", lwd=2, xlim=xlimits, ylim=ylimits, xlab="PC", cex=1.5, cex.axis=1.25, ylab=ylab, main=main, cex.lab=1.2)
#  grid(nx=NULL, ny=NA, lty=1)
  grid(nx=NA, ny=NULL, lty=2)
  for (k in 1:npcs) { lines(c(k,k), c(ylimits[1],100), col="grey") }
  abline(h=100, lwd=2)

  points(1:npcs, vars, pch=20, col="red")

  lines(1:npcs, varsLow, pch=20, type="b", cex=1.5, lwd=1, lty=2)
  points(1:npcs, varsLow, pch=20, col="red")
  lines(1:npcs, varsHigh, pch=20, type="b", cex=1.5, lwd=1, lty=2)
  points(1:npcs, varsHigh, pch=20, col="red")
}
if ((cumul==0)|(cumul==2)) {
  Fvars <- FRBres$eigval
  if (confmethod=="basic") {
    FvarsLow <- FRBres$eigval.CI.basic[,1]
    FvarsHigh <- FRBres$eigval.CI.basic[,2]
  } else {
    FvarsLow <- FRBres$eigval.CI.bca[,1]
    FvarsHigh <- FRBres$eigval.CI.bca[,2]
  }

  vars <- Fvars[1:npcs]
  varsLow <- FvarsLow[1:npcs]
  varsHigh <- FvarsHigh[1:npcs]

  xlimits <- c(0.5,npcs+0.5)
  if (cumul==2) {
    ylab=""
    main=paste("Variance (+ ",conf*100,"% ",confmethod," confidence limits)",sep="")
  }
  else {
    ylab="Variance"
    main=paste("PC variances (+ ",conf*100,"% ",confmethod," confidence limits",sep="")
  }
  plot(1:npcs, varsHigh, pch=20, type="n", xlim=xlimits, xlab="PC", cex=1.5, cex.axis=1.25, ylab=ylab, main=main, cex.lab=1.2)
  lines(1:npcs, vars, pch=20, type="b", cex=1.5, lwd=2)
#  grid(nx=NULL, ny=NA, lty=1)
  grid(nx=NA, ny=NULL, lty=2)
  for (k in 1:npcs) { abline(v=k, col="grey") }
  abline(h=0, lwd=2)

  points(1:npcs, vars, pch=20, col="red")

  lines(1:npcs, varsLow, pch=20, type="b", cex=1.5, lwd=1, lty=2)
  points(1:npcs, varsLow, pch=20, col="red")
  lines(1:npcs, varsHigh, pch=20, type="b", cex=1.5, lwd=1, lty=2)
  points(1:npcs, varsHigh, pch=20, col="red")
}

}

