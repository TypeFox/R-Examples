# this provides a modified version of the confint function from the bbmle package
# to fix bugs and increase functionality
# it is hoped to eventually roll these changes into bbmle

setMethod("confint", "profile.mymle",
function (object, parm, level = 0.95, trace=FALSE, ...)
{
  Pnames <- names(object@profile)
  if (missing(parm)) parm <- Pnames
  if (is.character(parm)) parm <- match(parm,Pnames)
  if (any(is.na(parm))) stop("parameters not found in profile")
  ## Calculate confidence intervals based on likelihood
  ## profiles
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(round(100 * a, 1), "%")
  ci <- array(NA, dim = c(length(parm), 2),
              dimnames = list(Pnames[parm], pct))
  cutoff <- qnorm(a)
  for (pm in parm) {
    pro <- object@profile[[Pnames[pm]]]
    pv <- pro[,"par.vals"]
    if (is.matrix(pv)) pv <- pv[,Pnames[pm]]
    if (any(diff(pro[,1])<0)) {
      warning(paste("non-monotonic profile (",
                    Pnames[pm],"): reverting from spline to linear approximation ",
              "(consider running 'profile' with manually reduced std.err)", sep=""))
      tt <- approx(pro[,1],pv,xout=cutoff)$y
    } else {
      sp <- spline(x = pv, y = pro[, 1])
      if (any(diff(sp$y)<0)) {
        warning(paste("non-monotonic spline fit to profile (",
                      Pnames[pm],"): reverting from spline to linear approximation",sep=""))
        tt <- approx(pro[,1],pv,xout=cutoff)$y
      } else {
        tt <- try(approx(sp$y, sp$x, xout = cutoff)$y,silent=TRUE)
        if (inherits(tt,"try-error")) tt <- rep(NA,2)
      }
    }
    ci[Pnames[pm], ] <- tt
  }
  drop(ci)
})

setMethod("confint", "mymle",
function (object, parm, level = 0.95, method,
          trace=FALSE,quietly=!interactive(),
          tol.newmin=0.001,...)
{
  if (missing(method)) method <- mymle.options("confint")
  ## changed coef() calls to object@coef -- really *don't* want fullcoef!
  Pnames <- names(object@coef)
  if (missing(parm))
    parm <- seq(along=Pnames)
  if (is.character(parm)) parm <- match(parm,Pnames)
  if (any(is.na(parm))) stop("parameters not found in model coefficients")
  if (method=="spline") {
    if (!quietly) message("Profiling...\n")
    newpars_found <- FALSE
    prof = try(profile(object,which=parm,tol.newmin=tol.newmin,...))
    if (inherits(prof,"try-error")) stop(paste("Problem with profiling:",prof))
    if (class(prof)=="mymle") newpars_found <- TRUE
    if (newpars_found) {
        ## profiling found a better fit
        message("returning better fit\n")
        return(prof)
    }
    return(confint(prof, parm, level, ...))
  } else {
    B0 <- object@coef
    pnames <- names(B0)
    if (missing(parm))
      parm <- seq(along=pnames)
    if (is.character(parm))
      parm <- match(parm, pnames, nomatch = 0)
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(round(100 * a, 1), "%")
    pct <- paste(round(100 * a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2),
                dimnames = list(pnames[parm], pct))
    std.err <- summary(object)@coef[, "Std. Error"]
    if (method=="uniroot") {
      chisqcutoff <- qchisq(level,1)
      call <- object@call
      if (!isTRUE(call$vecpar))
        call$start <- as.list(B0) ## added
      upper <- rep(unlist(eval.parent(call$upper)),length.out=length(pnames))
      lower <- rep(unlist(eval.parent(call$lower)),length.out=length(pnames))
      for (pm in parm) {
        critfun <- function(bi)
          {
            fix <- list(bi)
            names(fix) <- pnames[pm]
            call$fixed <- c(fix,eval(call$fixed))
            if (!is.null(upper) && length(upper)>1) call$upper <- upper[-pm]
            if (!is.null(lower) && length(lower)>1) call$lower <- lower[-pm]
            pfit <- try(eval(call), silent=TRUE)
            if(inherits(pfit, "try-error")) {
              warning(paste("Error encountered in profile (uniroot):",pfit))
              return(NA)
            }
            else {
              zz <- 2*pfit@min - 2*(-logLik(object))
              if (zz > -tol.newmin)
                zz <- max(zz, 0)
              else
                stop(sprintf("profiling has found a better solution (old deviance=%.2f, new deviance=%.2f), so original fit had not converged",2*pfit@min,2*(-logLik(object))))
              z <- zz - chisqcutoff
            }
            if (trace) cat(bi, z, "\n")
            z
          }
        stepfun <- function(step) {
          B0[pm] + sgn * step * std.err[pm]
        }
        invstepfun <- function(out) {
          (out - B0[pm])/(sgn * std.err[pm])
        }
        sgnvec=c(-1,1)
        for (i in 1:2) {
          sgn <- sgnvec[i]
          bnd <- if (sgn<0) {
            if (is.null(lower)) -Inf else lower[pm]
          } else {
            if (is.null(upper)) Inf else upper[pm]
          }
          c0 <- critfun(B0[pm])
          bi <- 
          ctry <- pmin(5,invstepfun(bnd))
          cdel <- -0.25
          c5 <- NA
          while (is.na(c5) && ctry>0 ) {
            c5 <- critfun(stepfun(ctry))
            if (is.na(c5)) {
              if (trace) cat("encountered NA, reducing ctry to",ctry+cdel,"\n")
              ctry <- ctry+cdel
            }
          }
          if (trace) cat(c0,c5,"\n")
          if (is.na(c0*c5) || c0*c5>0) {
            warning(paste("can't find confidence limits in",
                          c("negative","positive")[i],"direction"))
            curci <- NA
            ## FIXME: could try harder!
          } else {
            curci <- uniroot(critfun,c(stepfun(0),stepfun(ctry)))$root
          }
          ci[pnames[pm],i] <- curci
        }
      }
    } else if (method=="quad") {
      for (pm in parm) {
        ci[pnames[pm],] <- qnorm(a,B0[pm],std.err[pm])
      }
    } else stop("unknown method")
    return(drop(ci))
  }
})
