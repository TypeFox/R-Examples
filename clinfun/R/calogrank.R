calogrank <- function(ftime, fstatus, grp, cvt, strat=NULL) {
  requireNamespace("survival")
  call <- match.call()
  n0 <- length(ftime)
  cvt <- cbind(cvt)
  if (missing(strat)) strat <- rep(1,n0)
  if (length(fstatus) != n0 | length(grp) != n0 | nrow(cvt) != n0) stop("Data length mismatch")
  p <- ncol(cvt)
# if more than 3 covariates stop
  if (p > 3) stop("if more than 3 covariates use with principal components")

  ii <- !is.na(ftime) & !is.na(fstatus) & !is.na(grp) & !is.na(strat)
  for (i in 1:p) ii <- ii & !is.na(cvt[,i])
  ftime	<- ftime[ii]
  fstatus <- fstatus[ii]
  grp <- grp[ii]
  cvt <- cbind(cvt[ii,])
  strat <- strat[ii]
 
  n0 <- length(ftime)
  ng <- length(unique(grp))
  nstrat <- length(unique(strat))

  bb <- apply(cvt, 2, sd)/n0^0.26  

  rval <- survival::survdiff(Surv(ftime, fstatus) ~ grp + strata(strat))
  rval$call <- call

  ord <- order(strat, ftime, -fstatus)
  osts <- fstatus[ord]
  ogrp <- as.numeric(as.factor(grp))[ord]
  ocov <- as.matrix(cvt[ord,])
  strat <- strat[ord]
  ustat <- rep(0, ng)
  uvar <- matrix(0, ng, ng)
  uexpected <- rep(0, ng*nstrat)
  l <- 0
  for(i in unique(strat)) {
    n <- sum(strat==i)
    zz <- .Fortran("uclrst",
                   as.integer(n),
                   as.integer(ng),
                   as.integer(p),
                   as.double(osts[strat==i]),
                   as.integer(ogrp[strat==i]),
                   as.double(ocov[strat==i,]),
                   a0=double(n),
                   a1=double(n*ng),
                   xi=double(p),
                   xj=double(p),
                   Vii=double(ng),
                   Vij=double(ng),
                   Vji=double(ng),
                   Vjj=double(ng),
                   Vidot=double(ng*n),
                   Vdotj=double(ng*n),
                   Vijij=double(ng*ng),
                   igrp=double(ng),
                   jgrp=double(ng),
                   lrmn=double(ng),
                   lrvar=double(ng*ng),
                   as.double(bb))
    ustat <- ustat + zz$lrmn
    uvar <- uvar + matrix(zz$lrvar,ng,ng)
    uexpected[l*ng + (1:ng)] <- zz$jgrp
    l <- l+1
  }
  chi <- sum(solve(uvar[-1,-1], ustat[-1]) * ustat[-1])
  rval$uexp <- uexpected
  rval$uchisq <- chi
  rval$ustat <- ustat
  rval$uvar <- uvar
  class(rval) <- "calogrank"
  rval
}

print.calogrank <- function(x, digits = max(options()$digits - 4, 3), ...) {

  saveopt <-options(digits=digits)
  on.exit(options(saveopt))

  if (!inherits(x, 'calogrank')) stop("Object is not the result of calogrank")

  if (!is.null(cl<- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }

  omit <- x$na.action
  if (length(omit)) cat("n=", sum(x$n), ", ", naprint(omit),
					  ".\n\n", sep='')

  if (is.matrix(x$obs)){
    otmp <- apply(x$obs,1,sum)
    etmp <- apply(x$exp,1,sum)
    uetmp <- apply(x$uexp,1,sum)
  }
  else {
    otmp <- x$obs
    etmp <- x$exp
    uetmp <- x$uexp
  }
  df <- (sum(1*(etmp>0))) -1
  temp <- cbind(x$n, otmp, etmp, uetmp, ((otmp-etmp)^2)/ diag(x$var),
					  ((otmp-uetmp)^2)/ diag(x$uvar))
  dimnames(temp) <- list(names(x$n), c("N", "Observed", "Expected",
                                "Exptd-adj","(O-E)^2/V", "(O-aE)^2/aV"))
  print(temp)
  cat("\n Unconditional chisq=", format(round(x$chisq,2)),
      " on", df, "degrees of freedom, p=",
      format(signif(1-pchisq(x$chisq, df),digits)), "\n")
  
  cat("      Adjusted chisq=", format(round(x$uchisq,2)),
      " on", df, "degrees of freedom, p=",
      format(signif(1-pchisq(x$uchisq, df),digits)), "\n")
  
  invisible(x)
}
