"bugs.plot.inferences" <-
function (sims, display.parallel, ...){
  if (.Device=="windows" ||
      (.Device=="null device" && options("device")=="windows")){
    cex.names <- .7
    cex.axis <- .6
    cex.tiny <- .4
    cex.points <- .7
    standard.width <- 30
    max.width <- 40
    min.width <- .02
  }
  else {
    cex.names <- .7
    cex.axis <- .6
    cex.tiny <- .4
    cex.points <- .3
    standard.width <- 30
    max.width <- 40
    min.width <- .01
  }
  rootnames <- sims$root.short
  n.roots <- length(rootnames)
  sims.array <- sims$sims.array
  n.chains <- sims$n.chains
  dimension.short <- sims$dimension.short
  indexes.short <- sims$indexes.short
  long.short <- sims$long.short
  height <- .6
  par (mar=c(0,0,1,0))

  plot (c(0,1), c(-n.roots-.5,-.4),
        ann=FALSE, bty="n", xaxt="n", yaxt="n", type="n")

  W <- max(strwidth(rootnames, cex=cex.names))
  B <- (1-W)/3.8
  A <- 1-3.5*B
  if (display.parallel)
    text (A, -.4, "80% interval for each chain", adj=0, cex=cex.names)
  else
    text (A, -.4, "medians and 80% intervals", adj=0, cex=cex.names)
  num.height <- strheight (1:9, cex=cex.tiny)
  for (k in 1:n.roots){
    text (0, -k, rootnames[k], adj=0, cex=cex.names)
    J <- min (length(long.short[[k]]), max.width)
    if (k==1)
      index <- 1:J
    else
      index <- sum (unlist(lapply(long.short,length))[1:(k-1)]) + 1:J
    spacing <- 3.5/max(J,standard.width)
    med <- numeric(J)
    i80 <- matrix( , J, 2)
    med.chains <- matrix( , J, sims$n.chains)
    i80.chains <- array(NA, c(J, sims$n.chains, 2))
    for (j in 1:J){
      med[j] <- median(sims.array[,,index[j]])
      i80[j,] <- quantile(sims.array[,,index[j]], c(.1,.9))
      for (m in 1:n.chains){
        med.chains[j,m] <- quantile (sims.array[,m,index[j]], .5)
        i80.chains[j,m,] <- quantile (sims.array[,m,index[j]], c(.1,.9))
      }
    }
    rng <- range (i80, i80.chains)
    p.rng <- pretty(rng, n = 2)
    b <- height/(max(p.rng) - min(p.rng))
    a <- -(k+height/2) - b*p.rng[1]
    lines (A+c(0,0), -k+c(-height/2,height/2))
#
# plot a line at zero (if zero is in the range of the mini-plot)
#    
    if (min(p.rng)<0 & max(p.rng)>0)
      lines (A+B*spacing*c(0,J+1), rep (a,2), lwd=.5, col="gray")
    
    for (x in p.rng){
      text (A-B*.2, a+b*x, x, cex=cex.axis)
      lines (A+B*c(-.05,0), rep(a+b*x,2))
    }
    for (j in 1:J){
      if (display.parallel){
        for (m in 1:n.chains){
          interval <- a + b*i80.chains[j,m,]
          if (interval[2]-interval[1] < min.width)
            interval <- mean(interval) + c(-1,1)*min.width/2
          lines (A+B*spacing*rep(j+.6*(m-(n.chains+1)/2)/n.chains,2),
                 interval, lwd=.5, col=m+1)
        }
      }
      else {
        lines (A+B*spacing*rep(j,2), a + b*i80[j,], lwd=.5)
        for (m in 1:n.chains)
#        points (A+B*spacing*j, a + b*med[j], pch=20, cex=cex.points)
          points (A+B*spacing*j, a + b*med.chains[j,m], pch=20, cex=cex.points,
                  col=m+1)
      }
      dk <- dimension.short[k]
      if (dk>0){
        for (m in 1:dk){
          index0 <- indexes.short[[k]][[j]][m]
          if (j==1)
            text(A+B*spacing*j, -k-height/2-.05-num.height*(m-1), index0,
              cex=cex.tiny)
          else if (index0!=indexes.short[[k]][[j-1]][m] &
            (index0%%(floor(log10(index0)+1))==0))
              text(A+B*spacing*j, -k-height/2-.05-num.height*(m-1), index0,
                cex=cex.tiny)
        }
      }
    }
    if (J<length(long.short[[k]])) text (-.015, -k, "*",
                                         cex=cex.names, col="red")
  }
}
