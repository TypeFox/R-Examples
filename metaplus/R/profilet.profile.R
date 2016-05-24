##
## modified from profile method of bbmle package, which was modified from mle of the stats4 package
## 
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
## 
setGeneric("profilet.profile", function(fitted, which = 1:p, maxsteps = 100,
                                          alpha = 0.01, zmax = sqrt(qchisq(1 - alpha/2, p)),
                                          del = zmax/5, trace = FALSE, skiperrs=TRUE,
                                          std.err, tol.newmin = 0.001, debug=FALSE,
                                          prof.lower, prof.upper, skip.hessian=TRUE,
                                          try_harder=FALSE, ...)
                                          		{ standardGeneric("profilet.profile")})

setMethod("profilet.profile", "mymle",
          function (fitted, which = 1:p, maxsteps = 100,
                    alpha = 0.01, zmax = sqrt(qchisq(1 - alpha/2, p)),
                    del = zmax/5, trace = FALSE, skiperrs=TRUE,
                    std.err, tol.newmin = 0.001, debug=FALSE,
                    prof.lower, prof.upper, skip.hessian=TRUE,
                    try_harder=FALSE, ...) {
            ## fitted: mle2 object
            ## which: which parameters to profile (numeric or char)
            ## maxsteps: steps to take looking for zmax
            ## alpha: max alpha level
            ## zmax: max log-likelihood difference to search to
            ## del: stepsize
            ## trace:
            ## skiperrs:
             if (fitted@optimizer=="optimx") {
              fitted@call$method <- fitted@details$method.used
            }
            if (fitted@optimizer=="constrOptim")
              stop("profiling not yet working for constrOptim -- sorry")
            Pnames <- names(fitted@coef)
            p <- length(Pnames)
            if (is.character(which)) which <- match(which,Pnames)
            if (any(is.na(which)))
              stop("parameters not found in model coefficients")
            ## global flag for better fit found inside profile fit
            newpars_found <- FALSE
            if (debug) cat("i","bi","B0[i]","sgn","step","del","std.err[i]","\n")
            onestep <- function(step,bi) {
              if (missing(bi)) {
                bi <- B0[i] + sgn * step * del * std.err[i]
                if (debug) cat(i,bi,B0[i],sgn,step,del,std.err[i],"\n")
              } else if (debug) cat(bi,"\n")
              fix <- list(bi)
              names(fix) <- p.i
              if (is.null(call$fixed)) call$fixed <- fix
              else call$fixed <- c(eval(call$fixed),fix)
              save.start <- makestart.profilenorm.metaplus(call$data$yi,call$data$sei,mods=call$data$mods,fixed=call$fixed)
              maxll <- -Inf
              for (vinv in c(0.0,0.01,0.05,0.1,0.2,0.5,1)) { 
                start <- save.start
                if (length(start) >= 3) start <- c(start[1:2],vinv=vinv,start[3:length(start)])
                else start <- c(start,vinv=vinv)
                #browser()
                #start <- start[-i]
                call$start <- start
                #browser()
                if (skiperrs) {
                  pfit <- try(eval.parent(call, 2L), silent=TRUE)
                } else {
                  pfit <- eval.parent(call, 2L)
                }
                
                #print(profilet.fit)
                #browser()
                if (logLik(pfit) > maxll) {
                  maxfit <- pfit
                  maxll <- logLik(pfit)
                }
              }
              pfit <- maxfit              
              
              ok <- ! inherits(pfit,"try-error")
              if (debug && ok) cat(coef(pfit),-logLik(pfit),"\n")
              if(skiperrs && !ok) {
                warning(paste("Error encountered in profile:",pfit))
                return(NA)
              }
              else {
                ## pfit is current (profile) fit,
                ##   fitted is original fit
                ## pfit@min _should_ be > fitted@min
                ## thus zz below should be >0
                zz <- 2*(pfit@min - fitted@min)
                ri <- pv0
                ri[, names(pfit@coef)] <- pfit@coef
                ri[, p.i] <- bi
                ##cat(2*pfit@min,2*fitted@min,zz,
                ##   tol.newmin,zz<(-tol.newmin),"\n")
                if (!is.na(zz) && zz<0) {
                  if (zz > (-tol.newmin)) {
                    z <- 0
                    ## HACK for non-monotonic profiles? z <- -sgn*sqrt(abs(zz))
                  } else {
                    ## cat() instead of warning(); FIXME use message() instead???
                    message("Profiling has found a better solution,",
                            "so original fit had not converged:\n")
                    message(sprintf("(new deviance=%1.4g, old deviance=%1.4g, diff=%1.4g)",
                                    2*pfit@min,2*fitted@min,2*(pfit@min-fitted@min)),"\n")
                    #browser()
                    message("Returning better fit ...\n")
                    ## need to return parameters all the way up
                    ##   to top level
                    newpars_found <<- TRUE
                    ## return(pfit@fullcoef)
                    return(pfit) ## return full fit
                  }
                } else {
                  z <- sgn * sqrt(zz)
                }
                pvi <<- rbind(pvi, ri)
                zi <<- c(zi, z) ## nb GLOBAL set
              }
              if (trace) cat(bi, z, "\n")
              z
            } ## end onestep
            ## Profile the likelihood around its maximum
            ## Based on profile.glm in MASS
            #browser()
             summ <- summary(fitted)
            if (missing(std.err)) {
              std.err <- summ@coef[, "Std. Error"]
            } else {
              n <- length(summ@coef)
              if (length(std.err)<n)
                std.err <- rep(std.err,length.out=length(summ@coef))
              if (any(is.na(std.err)))
                std.err[is.na(std.err)] <- summ@coef[is.na(std.err)]
            }
            if (any(is.na(std.err))) {
              std.err[is.na(std.err)] <- sqrt(1/diag(fitted@details$hessian))[is.na(std.err)]
              if (any(is.na(std.err))) {  ## still bad
                #browser()
                stop("Hessian is ill-behaved or missing, ",
                     "can't find an initial estimate of std. error ",
                     "(consider specifying std.err in profile call)")
              }
              ## warn anyway ...
              warning("Non-positive-definite Hessian, ",
                      "attempting initial std err estimate from diagonals")
            }
            Pnames <- names(B0 <- fitted@coef)
            pv0 <- t(as.matrix(B0))
            p <- length(Pnames)
            prof <- vector("list", length = length(which))
            names(prof) <- Pnames[which]
            call <- fitted@call
            call$skip.hessian <- skip.hessian ## BMB: experimental
            call$minuslogl <- fitted@minuslogl
            ndeps <- eval.parent(call$control$ndeps)
            parscale <- eval.parent(call$control$parscale)
            nc <- length(fitted@coef)
            xf <- function(x) if (is.null(x)) NULL else rep(x,length.out=nc) ## expand to length
            upper <- xf(unlist(eval.parent(call$upper)))
            lower <- xf(unlist(eval.parent(call$lower)))
            if (all(upper==Inf & lower==-Inf)) {
              lower <- upper <- NULL
              ## kluge: lower/upper may have been set to +/- Inf
              ##   in previous rounds, 
              ##  but we don't want them in that case
            }
            if (!missing(prof.lower)) prof.lower <- xf(prof.lower)
            if (!missing(prof.upper)) prof.upper <- xf(prof.upper)
            ## cat("upper\n")
            ## print(upper)
            stop_msg <- list()
            for (i in which) {
              zi <- 0
              pvi <- pv0
              p.i <- Pnames[i]
              wfun <- function(txt) paste(txt," (",p.i,")",sep="")
              stop_msg[[i]] <- list(down="",up="")
              for (sgn in c(-1, 1)) {
                dir_ind <- (sgn+1)/2+1 ## (-1,1) -> (1,2)
                if (trace) {
                  cat("\nParameter:", p.i, c("down", "up")[dir_ind], "\n")
                  cat("par val","sqrt(dev diff)\n")
                }
                step <- 0
                z <- 0
                ## This logic was a bit frail in some cases with
                ## high parameter curvature. We should probably at least
                ## do something about cases where the mle2 call fails
                ## because the parameter gets stepped outside the domain.
                ## (We now have.)
                call$start <- as.list(B0)
                lastz <- 0
                valf <- function(b) {
                  (!is.null(b) && length(b)>1) ||
                    (length(b)==1 && i==1 && is.finite(b))
                }
                lbound <- if (!missing(prof.lower)) {
                  prof.lower[i]
                } else if (valf(lower))
                { lower[i]
                } else -Inf
                ubound <- if (!missing(prof.upper)) prof.upper[i] else if (valf(upper)) upper[i] else Inf
                stop_bound <- stop_na <- stop_cutoff <- stop_flat <- FALSE
                while ((step <- step + 1) < maxsteps &&
                         ## added is.na() test for try_harder case
                         ## FIXME: add unit test!
                         (is.na(z) || abs(z) < zmax)) {
                  curval <- B0[i] + sgn * step * del * std.err[i]
                  if ((sgn==-1 & curval<lbound) ||
                        (sgn==1 && curval>ubound)) {
                    stop_bound <- TRUE;
                    stop_msg[[i]][[dir_ind]] <- paste(stop_msg[[i]][[dir_ind]],wfun("hit bound"))
                    break
                  }
                  z <- onestep(step)
                  ## stop on flat spot, unless try_harder
                  if (step>1 && (identical(oldcurval,curval) || identical(oldz,z))) {
                    stop_flat <- TRUE
                    stop_msg[[i]][[dir_ind]] <- paste(stop_msg[[i]][[dir_ind]],wfun("hit flat spot"),
                                                      sep=";")
                    if (!try_harder) break
                  }
                  oldcurval <- curval
                  oldz <- z
                  if (newpars_found) return(z)
                  if(is.na(z)) {
                    stop_na <- TRUE
                    stop_msg[[i]][[dir_ind]] <- paste(stop_msg[[i]][[dir_ind]],wfun("hit NA"),sep=";")
                    if (!try_harder) break
                  }
                  lastz <- z
                  if (newpars_found) return(z)
                }
                stop_cutoff <- (!is.na(z) && abs(z)>=zmax)
                stop_maxstep <- (step==maxsteps)
                if (stop_maxstep) stop_msg[[i]][[dir_ind]] <- paste(stop_msg[[i]][[dir_ind]],wfun("max steps"),sep=";")
                if (debug) {
                  if (stop_na) message(wfun("encountered NA"),"\n")
                  if (stop_cutoff) message(wfun("above cutoff"),"\n")
                }
                if (stop_flat) {
                  warning(wfun("stepsize effectively zero/flat profile"))
                } else {
                  if (stop_maxstep) warning(wfun("hit maximum number of steps"))
                  if(!stop_cutoff) {
                    if (debug) cat(wfun("haven't got to zmax yet, trying harder"),"\n")
                    stop_msg[[i]][[dir_ind]] <- paste(stop_msg[[i]][[dir_ind]],wfun("past cutoff"),sep=";")
                    ## now let's try a bit harder if we came up short
                    for(dstep in c(0.2, 0.4, 0.6, 0.8, 0.9)) {
                      curval <- B0[i] + sgn * (step-1+dstep) * del * std.err[i]
                      if ((sgn==-1 & curval<lbound) ||
                            (sgn==1 && curval>ubound)) break
                      z <- onestep(step - 1 + dstep)
                      if (newpars_found) return(z)
                      if(is.na(z) || abs(z) > zmax) break
                      lastz <- z
                      if (newpars_found) return(z)
                    }
                    if (!stop_cutoff && stop_bound) {
                      if (debug) cat(wfun("bounded and didn't make it, try at boundary"),"\n")
                      ## bounded and didn't make it, try at boundary
                      if (sgn==-1 && B0[i]>lbound) z <- onestep(bi=lbound)
                      if (sgn==1  && B0[i]<ubound) z <- onestep(bi=ubound)
                      if (newpars_found) return(z)
                    }
                  } else if (length(zi) < 5) { # try smaller steps
                    if (debug) cat(wfun("try smaller steps"),"\n")
                    stop_msg[[i]][[dir_ind]] <- paste(stop_msg[[i]][[dir_ind]],wfun("took more steps"),sep=";")
                    mxstep <- step - 1
                    step <- 0.5
                    while ((step <- step + 1) < mxstep) {
                      z <- onestep(step)
                    }
                  } ## smaller steps
                } ## !zero stepsize
              } ## step in both directions
              si <- order(pvi[, i])
              prof[[p.i]] <- data.frame(z = zi[si])
              prof[[p.i]]$par.vals <- pvi[si,, drop=FALSE]
            } ## for i in which
            newprof <- new("profile.mymle", profile = prof, summary = summ)
            attr(newprof,"stop_msg") <- stop_msg
            newprof
          })
