"tailordCOP" <-
function(cop=NULL, para=NULL, tol=1e-6, plot=FALSE, verbose=FALSE,  ...) {
  resolution <- abs(log10(1E-6)) - 1
  kapu.tmp <- kapu <- NA; u <- 0
  KAPU1 <- KAPU2 <- vector(mode="numeric")
  if(verbose) message("Upper-Tail Order")
  t <- 0.5
  while(t > .Machine$double.eps) {
    opts <- options(warn=-1)
    kapu.tmp <- log(surCOP(t,t, cop=cop, para=para, ...))/log(t)
    options(opts)
    if(is.nan(kapu.tmp) | ! is.finite(kapu.tmp)) break
    u <- u + 1
    if(verbose) message("  iter.",u,"  t=",t,"  kappaU=",kapu)
    KAPU1[u] <- t; KAPU2[u] <- kapu
    if(! is.na(kapu) & abs(kapu - kapu.tmp) < tol) break
    t <- t - t/2
    kapu <- kapu.tmp
  }
  if(u == 0) {
    if(! verbose) warning("  Apparent infinity on Upper-Tail Order, setting t=1/10000 and recomputing")
    t <- 0.0001; u <- 1
    kapu <- log(surCOP(t,t, cop=cop, para=para, ...))/log(t)
    KAPU1[u] <- t; KAPU2[u] <- kapu
    if(verbose) message("  iter.",u,"  t=",t,"  kappaU=",kapu)
  }

  if(verbose) message("\nLower-Tail Order")
  kapl.tmp <- kapl <- NA; l <- 0
  KAPL1 <- KAPL2 <- vector(mode="numeric")
  t <- 0.5
  while(t > .Machine$double.eps) {
    opts <- options(warn=-1)
    kapl.tmp <- log(cop(t,t, para=para, ...))/log(t)
    options(opts)
    if(is.nan(kapl.tmp) | ! is.finite(kapl.tmp)) break
    l <- l + 1
    if(verbose) message("  iter.",l,"  t=",t,"  kappaL=",kapl)
    KAPL1[l] <- t; KAPL2[l] <- kapl
    if(! is.na(kapl) & abs(kapl - kapl.tmp) < tol) break
    t <- t - t/2
    kapl <- kapl.tmp
  }
  if(l == 0) {
    if(! verbose) warning("  Apparent infinity on Lower-Tail Order, setting t=1/10000 and recomputing")
    t <- 0.0001; l <- 1
    kapl <- log(COP(t,t, cop=cop, para=para, ...))/log(t)
    KAPL1[l] <- t; KAPL2[l] <- kapl
    if(verbose) message("  iter.",l,"  t=",t,"  kappaL=",kapl)
  }
  if(verbose) message("---done.")
  if(plot) {
    xmin <- qnorm(min(1-KAPU1,KAPL1, na.rm=TRUE))
    xmax <- qnorm(max(1-KAPU1,KAPL1, na.rm=TRUE))
    ymin <- min(KAPU2,KAPL2, na.rm=TRUE)
    ymax <- max(KAPU2,KAPL2, na.rm=TRUE)
    plot(qnorm(1-KAPU1), KAPU2, type="l", col=4,
           xlim=c(xmin, xmax), ylim=c(ymin, ymax),
           xlab="STANDARD NORMAL VARIATES",
           ylab="UPPER- AND LOWER-TAIL ORDER")
    lines(qnorm(KAPL1), KAPL2, col=2)
    points(qnorm(1-KAPU1[u]), KAPU2[u], cex=0.7, pch=16, col=4)
    points(qnorm(1-KAPU1[u]), KAPU2[u], pch="U", cex=1.3)
    points(qnorm(KAPL1[l]), KAPL2[l], cex=0.7, pch=16, col=2)
    points(qnorm(KAPL1[l]), KAPL2[l], pch="L", cex=1.3)
  }
  return(list(kappaL = round(KAPL2[l], digits=resolution),
              kappaU = round(KAPU2[u], digits=resolution),
              source = "tailordCOP"))
}

