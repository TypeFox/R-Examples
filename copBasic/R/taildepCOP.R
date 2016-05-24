"taildepCOP" <-
function(cop=NULL, para=NULL, tol=1e-6, divisor=2, plot=FALSE, ylim=NULL, verbose=FALSE, ...) {
  resolution <- abs(log10(1E-6)) - 1
  lamu.tmp <- lamu <- u <- 0
  LAMU1 <- LAMU2 <- vector(mode="numeric")
  if(verbose) message("Upper-Tail Dependency")
  t <- 0.5
  while(t < 1 - .Machine$double.eps) {
    lamu.tmp <- 2 - (1-cop(t,t, para=para, ...))/(1-t)
    if(is.nan(lamu.tmp)) break
    u <- u + 1
    if(verbose) message("  iter.",u,"  t=",t,"  lambdaU=",lamu)
    LAMU1[u] <- t; LAMU2[u] <- lamu
    if(abs(lamu - lamu.tmp) < tol) break
    t <- t + (1-t)/divisor
    lamu <- lamu.tmp
  }
  if(verbose) message("\nLower-Tail Dependency")
  laml.tmp <- laml <- l <- 0
  LAML1 <- LAML2 <- vector(mode="numeric")
  t <- 0.5
  while(t > .Machine$double.eps) {
    laml.tmp <- cop(t,t, para=para, ...)/t
    if(is.nan(laml.tmp)) break
    l <- l + 1
    if(verbose) message("  iter.",l,"  t=",t,"  lambdaL=",laml)
    LAML1[l] <- t; LAML2[l] <- laml
    if(abs(laml - laml.tmp) < tol) break
    t <- t - t/divisor
    laml <- laml.tmp
  }
  if(verbose) message("---done.")
  if(plot) {
    xmin <- qnorm(min(LAMU1,LAML1))
    xmax <- qnorm(max(LAMU1,LAML1))
    ymin <- min(LAMU2,LAML2)
    ymax <- max(LAMU2,LAML2)
    if(is.null(ylim)) ylim <- c(ymin, ymax)
    plot(qnorm(LAMU1), LAMU2, type="l", col=4,
           xlim=c(xmin, xmax), ylim=ylim,
           xlab="STANDARD NORMAL VARIATES",
           ylab="UPPER- AND LOWER-TAIL DEPENDENCY PARAMETER")
    lines(qnorm(LAML1), LAML2, col=2)
    points(qnorm(LAMU1[u]), LAMU2[u], cex=0.7, pch=16, col=4)
    points(qnorm(LAMU1[u]), LAMU2[u], pch="U", cex=1.3)
    points(qnorm(LAML1[l]), LAML2[l], cex=0.7, pch=16, col=2)
    points(qnorm(LAML1[l]), LAML2[l], pch="L", cex=1.3)
  }
  return(list(lambdaL = round(LAML2[l], digits=resolution),
              lambdaU = round(LAMU2[u], digits=resolution),
              source = "taildepCOP"))
}

