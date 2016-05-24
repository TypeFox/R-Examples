"MM.E.gauss" <-
function(y, cov=FALSE, isigma=1, onlyS=FALSE, test=FALSE, level=0.90, control, ...){
  call <- match.call(); n   <- length(y)
  old <- comval();      tlo <- control$tlo; mxf <- control$mxf
  mxs <- control$mxs;   ntm <- control$ntm; tls <- control$tls  
  k0  <- control$k0;    k1  <- control$k1 ; h   <- control$h
  Beta<- control$Beta0 
  zc  <- list(...)
  if (length(zc) != 0) {
    if (!is.null(zc$tlo))    tlo <- zc$tlo
    if (!is.null(zc$mxf))    mxf <- zc$mxf
    if (!is.null(zc$mxs))    mxs <- zc$mxs
    if (!is.null(zc$ntm))    ntm <- zc$ntm
    if (!is.null(zc$tls))    tls <- zc$tls
    if (!is.null(zc$k0 ))    k0  <- zc$k0
    if (!is.null(zc$k1 ))    k1  <- zc$k1
    if (!is.null(zc$h  ))    h   <- zc$h
    if (!is.null(zc$Beta0)) Beta <- zc$Beta0}

# 1. Initial estimates (l0,s0) of location and scale 

  l0   <- median(y); s0 <- median(abs(y-l0))/.6745
  dfcomn2(ipsi=4, xk=k0, beta=Beta)
  y1   <- min(y); yn <- max(y)
  yseq <- matrix(seq(y1, yn, length=h+1), ncol=1)
  z    <- apply(yseq,1,s.estim,s0,y,tlo,mxs)
  z1   <- z[1,]
  zind <- seq(h+1)[z1==min(z1)]
  imin <- if (length(zind) > 1) sample(zind,1) else zind
  s0   <- z[1,imin]; l0 <- z[2,imin]
  if(s0 < tls) {cat(paste("Improved scale less than tls=", tls, "\n"))
                return(MM.E.null())}

# 2. Refinement

  zr <- lywalg(y, lambda=l0, sigmai=s0, tol=tlo, gam=1, 
               isigma=-1, maxit=mxf, maxis=1, nitmon=ntm)
  s1 <- zr$sigmaf; l1 <- zr$lambda
  if (zr$nit==mxf) cat("Maxit for refinement reached \n")
  if(s1 < s0) {l0 <- l1; s0 <- s1; r0 <- y - l0 } else 
  cat("Refinement of initial estimates did not decrease scale \n")
  if (onlyS) {ans <- list(lambda=l0,sigma=s0,tbias=NA,pchi=NA,V.lambda=NA,
             T.lambda=NA,T.sigma=NA,nit.ref=zr$nit,call=call); return(ans)}

# 3. Final estimates

  dfcomn2(ipsi=4, xk=k1)
  Qmm0 <- sum(Rho(r0/s0))/(n-1)
  zf <- lywalg(y, lambda=l0, sigmai=s0, tol=tlo, gam=1, 
               isigma=0, maxit=mxf, maxis=1, nitmon=ntm)
  l1 <- zf$lambda
  ans$nit.l1 <- zf$nit
  if(zf$nit==mxf) cat("Maxit for final estimate reached \n")
  r1   <- y-l1
  Qmm1 <- sum(Rho(r1/s0))/(n-1)
  if (Qmm1 > Qmm0+tlo) {cat("Step 3 does not decrease objective function: ")
                        cat("Qmm0=", Qmm0, " Qmm1=", Qmm1, "\n")}

# 4. Scale estimate for test

    dfcomn2(ipsi = 4, xk = k0, beta = Beta)
    zs <- rysigm(r1,wgt=y,sigmai=2*s0,np=1,itype=1,isigma=1,tol=tlo,maxis=mxs)
    s1 <- zs$sigmaf
    if(s1 < tls) {cat(paste("Final scale less than tls = ", tls, "\n")); return(MM.E.null())}
    ans$nit.sig <- zs$nit
    if(zs$nit == mxs) cat("Maxit for final sigma reached \n")

# 5. Test for bias

  tbias <- NA;  pchi  <- NA
  if (test) {
    dfcomn2(ipsi = 4, xk = k1)
    tmp <- r0/s0
    sc0 <- Psi(tmp)
    s1p <- sum(Psp(tmp))/n
    dfcomn2(ipsi = 4, xk = k0)
    s0p <- sum(Psp(tmp))/n
    sc1 <- Psi(tmp)
    tmp <- tmp * Psi(tmp)
    s0r <- (sum(tmp) * s0)/n
    if(s0r < tls | s0p < tls | s1p < tls) {
    cat(paste("Denominator smaller than tls=",tls," in test for bias\n")); return(MM.E.null())}
    tmp <- sc0/s1p - sc1/s0p
    d2  <- sum(tmp * tmp)/n
    if(d2 < tls) {
    cat(paste("Denominator smaller than tls=", tls," in test for bias\n")); return(MM.E.null())}
    tbias <- (2*n*(s1 - s0)*s0r)/(s0p*d2*s0*s0)
    qchi  <- qchisq(level, 1)
    pchi  <- pchisq(tbias, 1)
    if(tbias > qchi) cat(paste("Significant test at level ",100*level, "%\n"))}
    ans <- list(lambda=l1, sigma=s0, tbias=tbias, pchi=pchi)
# 6. Covariances

  if(cov) {z <- AV.MM.gauss(k0=k0,k1=k1, sigma=s0)
           ans$V.lambda <- z$V.lambda; ans$V.sigma <- z$V.sigma }
  ans$lambda <- l1; ans$T.lambda <- l0
  if (test) {
    if (tbias > qchi) {
     ans$lambda <- l0; ans$T.lambda <- l1; ans$V.lambda <- NA
     cat("\nTest for bias = ", format(tbias), " with 1",
         " degree of freedom. \nProbability( chi-squared > ",
           format(tbias), ") =", format(1-pchi), "\n")
     cat("\nWarning: the bias is high; inference based on final estimate is not\n")
     cat("         recommended; use initial estimate as exploratory tools.\n")}
  }
  ans$T.sigma <- s1 

# 7. Compute Qn scale estimate

  qn <- NA
  if (isigma==2) {qn <- Qn(y)$scale; ans$sigma <- qn
                  if (cov) ans$V.sigma <- qn*qn*0.6089}
  dfcomn2(ipsi = old$ipsi, xk = old$xk, beta = old$bta)
  ans$call <- call; ans}

