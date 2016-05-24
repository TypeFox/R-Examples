# this provides a modified version of the profile function from the bbmle package
# to fix bugs and increase functionality
# it is hoped to eventually roll these changes into bbmle

## FIXME: abstract to general-purpose code?  (i.e. replace 'fitted' by
#    objective function, parameter vector, optimizer, method, control settings,
##   min val, standard error/Hessian, ...
##
## allow starting values to be set by "mle" (always use mle), "prevfit"
##  (default?), and "extrap" (linear extrapolation from previous two fits)
## 

setMethod("profile", "mymle",
          function (fitted, which = 1:p, maxsteps = 100,
                    alpha = 0.01, zmax = sqrt(qchisq(1 - alpha/2, p)),
                    del = zmax/5, trace = FALSE, skiperrs=TRUE,
                    std.err, tol.newmin = 0.001, debug=FALSE,
                    prof.lower, prof.upper, skip.hessian=TRUE,
                    try_harder=FALSE, start.method, ...) {
              ## fitted: mle2 object
              ## which: which parameters to profile (numeric or char)
              ## maxsteps: steps to take looking for zmax
              ## alpha: max alpha level
              ## zmax: max log-likelihood difference to search to
              ## del: stepsize
              ## trace:
              ## skiperrs:
            if (missing(start.method)) start.method <- "prevfit"
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
                #browser()
                #print(as.vector(call$start))
                if (skiperrs) {
                    pfit <- try(eval.parent(call, 2L), silent=TRUE)
                } else {
                    pfit <- eval.parent(call, 2L)
                }
                ok <- ! inherits(pfit,"try-error")
                #print(pfit)
                #browser()
                #print(pfit@details$convergence)
                #print(pfit@details$message)
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
                if ((step>=1) & (start.method=="prevfit")) call$start <<- as.list(pfit@coef)                     
                z
            } ## end onestep
            ## Profile the likelihood around its maximum
            ## Based on profile.glm in MASS
            summ <- summary(fitted)
            #browser()
            if (missing(std.err)) {
                std.err <- summ@coef[, "Std. Error"]
            } else {
                #browser()
                n <- dim(summ@coef)[1]
                if (length(std.err)!=n) stop("length standard errors not equal to coefficients length")
# not certain what this was supposed to do - better to stop
#                 if (length(std.err)<n)
#                   std.err <- rep(std.err,length.out=length(summ@coef))
#                 if (any(is.na(std.err)))
#                   std.err[is.na(std.err)] <- summ@coef[is.na(std.err)]
            }
            if (any(is.na(std.err))) {
              std.err[is.na(std.err)] <- sqrt(1/diag(fitted@details$hessian))[is.na(std.err)]
              if (any(is.na(std.err))) {  ## still bad
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
              ## omit values from control vectors:
              ##   is this necessary/correct?
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


# setMethod("plot", signature(x="profile.mymle", y="missing"),
#           function (x, levels, which=1:p, conf = c(99, 95, 90, 80, 50)/100,
#           plot.confstr = TRUE, confstr = NULL, absVal = TRUE, add = FALSE,
#           col.minval="green", lty.minval=2,
#           col.conf="magenta", lty.conf=2,
#           col.prof="blue", lty.prof=1,
#           xlabs=nm, ylab="z",
#           onepage=TRUE,
#           ask=((prod(par("mfcol")) < length(which)) && dev.interactive() &&
#                !onepage),
#           show.points=FALSE,
#           main, xlim, ylim, ...)
# {
#     ## Plot profiled likelihood
#     ## Based on profile.nls (package stats)
#     obj <- x@profile
#     nm <- names(obj)
#     p <- length(nm)
#     ## need to save these for i>1 below
#     no.xlim <- missing(xlim)
#     no.ylim <- missing(ylim)    
#     if (is.character(which)) which <- match(which,nm)
#     ask_orig <- par(ask=ask)
#     op <- list(ask=ask_orig)
#     if (onepage) {
#         nplots <- length(which)
#         ## Q: should we reset par(mfrow), or par(mfg), anyway?
#         if (prod(par("mfcol")) < nplots) {
#             rows <- ceiling(round(sqrt(nplots)))
#             columns <- ceiling(nplots/rows)
#             mfrow_orig <- par(mfrow=c(rows,columns))
#             op <- c(op,mfrow_orig)
#           }
#       }
#     on.exit(par(op))
#     confstr <- NULL
#     if (missing(levels)) {
#         levels <- sqrt(qchisq(pmax(0, pmin(1, conf)), 1))
#         confstr <- paste(format(100 * conf), "%", sep = "")
#     }
#     if (any(levels <= 0)) {
#         levels <- levels[levels > 0]
#         warning("levels truncated to positive values only")
#     }
#     if (is.null(confstr)) {
#         confstr <- paste(format(100 * pchisq(levels^2, 1)), "%", sep = "")
#     }
#     mlev <- max(levels) * 1.05
#     ##    opar <- par(mar = c(5, 4, 1, 1) + 0.1)
#     if (!missing(xlabs) && length(which)<length(nm)) {
#       xl2 = nm
#       xl2[which] <- xlabs
#       xlabs <- xl2
#     }
#     if (missing(main)) 
#       main <- paste("Likelihood profile:",nm)
#     main <- rep(main,length=length(nm))
#     for (i in seq(along = nm)[which]) {
#         ## <FIXME> This does not need to be monotonic
#         ## cat("**",i,obj[[i]]$par.vals[,i],obj[[i]]$z,"\n")
#         ## FIXME: reconcile this with confint!
#         yvals <- obj[[i]]$par.vals[,nm[i],drop=FALSE]
#         avals <- data.frame(x=unname(yvals), y=obj[[i]]$z)
#         if (!all(diff(obj[[i]]$z)>0)) {
#           warning("non-monotonic profile: reverting to linear interpolation.  Consider setting std.err manually")
#           predback <- approxfun(obj[[i]]$z,yvals)
#         } else {
#           sp <- splines::interpSpline(yvals, obj[[i]]$z,
#                                       na.action=na.omit)
#           avals <- rbind(avals,as.data.frame(predict(sp)))
#           avals <- avals[order(avals$x),]
#           bsp <- try(splines::backSpline(sp),silent=TRUE)
#           bsp.OK <- (class(bsp)[1]!="try-error")
#           if (bsp.OK) {
#             predback <- function(y) { predict(bsp,y)$y }
#           } else { ## backspline failed
#             warning("backspline failed: using uniroot(), confidence limits may be unreliable")
#             ## what do we do?
#             ## attempt to use uniroot
#             predback <- function(y) {
#               pfun0 <- function(z1) {
#                 t1 <- try(uniroot(function(z) {
#                   predict(sp,z)$y-z1
#                 }, range(obj[[i]]$par.vals[,nm[i]])),silent=TRUE)
#                 if (class(t1)[1]=="try-error") NA else t1$root
#                 }
#               sapply(y,pfun0)
#             }
#           }
#         }
#         ## </FIXME>
#         if (no.xlim) xlim <- predback(c(-mlev, mlev))
#         xvals <- obj[[i]]$par.vals[,nm[i]]
#         if (is.na(xlim[1]))
#           xlim[1] <- min(xvals)
#         if (is.na(xlim[2]))
#           xlim[2] <- max(xvals)
#         if (absVal) {
#             if (!add) {
#                 if (no.ylim) ylim <- c(0,mlev)
#                 plot(abs(obj[[i]]$z) ~ xvals, 
#                      xlab = xlabs[i],
#                      ylab = if (missing(ylab)) expression(abs(z)) else ylab,
#                      xlim = xlim, ylim = ylim,
#                      type = "n", main=main[i], ...)
#             }
#             avals$y <- abs(avals$y)
#             lines(avals, col = col.prof, lty=lty.prof)
#             if (show.points) points(yvals,abs(obj[[i]]$z))
#         } else { ## not absVal
#             if (!add) {
#                 if (no.ylim) ylim <- c(-mlev,mlev)
#                 plot(obj[[i]]$z ~ xvals,  xlab = xlabs[i],
#                      ylim = ylim, xlim = xlim,
#                      ylab = if (missing(ylab)) expression(z) else ylab,
#                      type = "n", main=main[i], ...)
#             }
#             lines(avals, col = col.prof, lty=lty.prof)
#             if (show.points) points(yvals,obj[[i]]$z)
#         }
#         x0 <- predback(0)
#         abline(v = x0, h=0, col = col.minval, lty = lty.minval)
#         for (j in 1:length(levels)) {
#             lev <- levels[j]
#             confstr.lev <- confstr[j]
#             ## Note: predict may return NA if we didn't profile
#             ## far enough in either direction. That's OK for the
#             ## "h" part of the plot, but the horizontal line made
#             ## with "l" disappears.
#             pred <- predback(c(-lev, lev))
#             ## horizontal
#             if (absVal) levs=rep(lev,2) else levs=c(-lev,lev)
#             lines(pred, levs, type = "h", col = col.conf, lty = 2)
#             ## vertical
#             pred <- ifelse(is.na(pred), xlim, pred)
#             if (absVal) {
#                 lines(pred, rep(lev, 2), type = "l", col = col.conf, lty = lty.conf)
#             } else {
#                 lines(c(x0,pred[2]), rep(lev, 2), type = "l", col = col.conf, lty = lty.conf)
#                 lines(c(pred[1],x0), rep(-lev, 2), type = "l", col = col.conf, lty = lty.conf)
#             }
#             if (plot.confstr) {
#                 text(labels=confstr.lev,x=x0,y=lev,col=col.conf)
#             }
#         } ## loop over levels
#     } ## loop over variables
#     ## par(opar)
#   })
