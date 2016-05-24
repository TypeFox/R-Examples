## file cond/R/cond.R, v 1.2-3 2014-06-27
##
##  Copyright (C) 2000-2014 Alessandra R. Brazzale 
##
##  This file is part of the "cond" package for R.  This program is 
##  free software; you can redistribute it and/or modify it under the 
##  terms of the GNU General Public License as published by the Free 
##  Software Foundation; either version 2 of the License, or (at your 
##  option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
##  MA 02111-1307 USA or look up the web page 
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Alessandra R. Brazzale, Department of Statistics, University of 
##  Padova, Via C. Battisti 241/243, 35121 Padova (PD), Italy.
##  Email: alessandra.brazzale@unipd.it  
##  Web: http://www.stat.unipd/~brazzale

cond <- function(object, offset, ...)  UseMethod("cond")

cond.glm <- function(object, offset, formula = NULL, family = NULL,   
                     data = sys.frame(sys.parent()), pts = 20, 
                     n = max(100, 2*pts), tms = 0.6, from = NULL, 
                     to = NULL, control = glm.control(...),  
                     trace = FALSE, ...)
{
  m <- match.call()
  if(missing(offset))  
    stop("Argument \"offset\" is missing, with not default")
  if(missing(object))
  {
    if(is.null(formula) || is.null(family))
       stop("Model is missing, with no default")
    else 
    {
      new.call <- m
      new.call$offset <- new.call$pts <- new.call$n <- NULL
      new.call$tms <- new.call$from <- new.call$to <- NULL
      new.call$trace <- NULL
      new.call[[1]] <- as.name("glm")
      object <- eval(new.call, parent.frame()) 	
    }
  }
  else
  {
    oc <- attr(object, "class")
    if(is.null(oc) || oc[1] != "glm")
      stop("Invalid argument: not a \"glm\" object")
  }
  formula <- formula(object$terms)
  data <- object$data
  mf <- model.frame.default(object, data=data)
  Terms <- attr(mf, "terms")
  is.empty <- is.empty.model(Terms)
  if(is.empty)
    stop("Invalid model: all parameters are fixed")
  is.scalar <- (dim(model.matrix(object))[2] == 1)
  offsetName <- if(!is.character(m$offset)) deparse(m$offset,
                                                    width.cutoff = 500)
                else m$offset                
  .offsetName <- offsetName   
  if((offsetName == "1") && (attr(Terms, "intercept") == 0))
    stop(paste("Invalid argument for \"offset\": intercept not included in original model"))
  if(!match(offsetName, 
            c("1", if(!is.empty)
                     dimnames(model.matrix(object))[[2]]),  
            nomatch=FALSE))
    stop(paste("Invalid argument for \"offset\": variable not included in original model"))
  switch(offsetName, 
         "1"= { 
                nobs <- dim(as.matrix(model.response(mf)))[1]
                .offset <- rep(1, nobs) 
              }, 
              { 
                if(is.call(m$offset))
                {
                  offsetFactors <- attr(Terms, "factors")
                  if(sum(offsetFactors[,offsetName]) > 1)
                    stop("Invalid argument for \"offset\": parameter of interest cannot be an interaction")
                  else
                  .offset <- eval(m$offset, envir=data)
                }
                else
                {  
                  if(!match(offsetName, tl <- 
                            attr(Terms, "term.labels"), nomatch=FALSE))
                    .offsetName <- tl[as.logical(pmatch(tl, 
                                        offsetName, nomatch=FALSE))]
                  .offset <- do.call("eval", 
                                     list(expr=parse(text=.offsetName), 
                                          envir=data))
                  if(is.factor(.offset))
                  { 
                    offsetContrasts <- contrasts(.offset)
                    if(any(is.character(offsetContrasts)))
                      stop("Invalid argument for \"offset\": parameter of interest is not scalar")
                    if(any(dim(offsetContrasts) != c(2,1)))
                      stop("Invalid argument for \"offset\": parameter of interest is not scalar")
                    off.tmp <- numeric(length=length(.offset))
                    offsetLevels <- dimnames(offsetContrasts)[[1]]
                    off.tmp[.offset==offsetLevels[1]] <- 
                      offsetContrasts[1,]
                    off.tmp[.offset==offsetLevels[2]] <- 
                      offsetContrasts[2,]
                    if(any(is.na(.offset)))
                      off.tmp[is.na(.offset)] <- NA
                    .offset <- off.tmp
                  }
                }
              })              
#  assign(".offset", .offset, envir=sys.frame())
#  assign(".offset", .offset, pos=1)		## 20.05.13
#  on.exit(remove(".offset", pos=1), add=TRUE)
  if(!missing(pts) && (pts < 0))
    stop("Invalid argument: negative values not allowed for \"pts\"")
  else if(pts < 10)
          stop("Invalid argument: \"pts\" too small (< 10)")
  if(!missing(n) && (n < 0))
    stop("Invalid argument: negative values not allowed for \"n\"")
  else if(n < 50)
         stop("Invalid argument \"n\" too small (< 50)")
  if(!missing(tms) && (tms < 0))  
    stop("Invalid argument: negative values not allowed for \"tms\"")
  else if(tms > 1)
         warning("\"tms\" may be too large ( > 1)")
  summary.obj <- summary(object)
  glmMLE <- switch(offsetName,
                   "1" = summary.obj$coef["(Intercept)", "Estimate"],
                         summary.obj$coef[offsetName, "Estimate"])
  glmSE <- switch(offsetName,
                   "1" = summary.obj$coef["(Intercept)", "Std. Error"],
                         summary.obj$coef[offsetName, "Std. Error"])
  glmDeviance <- deviance(object)
#  glmDet <- 1/det(summary.obj$cov.unscaled)
  glmDet <- - determinant(summary.obj$cov.unscaled)$mod
  if(is.null(from))  
    from <- glmMLE - 3.5 * glmSE
  if(is.null(to))  
    to <- glmMLE + 3.5 * glmSE
  if(from > to)
    stop("Invalid range: \"from\" < \"to\"")
  pts <- pts + (pts%%2)
  offsetCoef <- seq(from, to, length=pts)
  lengthOC <- length(offsetCoef)
  glmCDev <- vector("numeric", lengthOC)
  if(!is.scalar)
    glmCDet <- glmCDev
#  if( (length(attr(Terms, "term.labels")) == 0) &&
#      !is.null(attr(Terms, "offset")) )
#  {
#    ..offset <- modOff <- model.offset(model.frame(object))
##    assign("..offset", ..offset, envir=sys.frame())
##    assign("..offset", ..offset, pos=1)		## 20.05.13
##    on.exit(remove("..offset", pos=1), add=TRUE)
#  }
  for(i in seq(along = offsetCoef))
  {
    if(trace) 
      cat(as.name(paste("\niteration", i,
                        switch(offsetName, 
                               "1"     = ": (Intercept)", 
                                paste(": coefficient of", offsetName)),
                        "=", signif(offsetCoef[i], digits=4))))
    .object <- object$call
    .object$formula <- formula
    switch(offsetName, 
           "1"     = { nf <- as.formula(".~. -1") },
                     { nf <- as.formula(paste(".~.-", .offsetName)) } )               
    nf <- update(.object$formula, nf)   
    .object$formula <- nf
    if(!is.null(.object$offset))
      .object$offset <- object$offset + offsetCoef[i]*.offset
    else
      .object$offset <- offsetCoef[i]*.offset  
    environment(.object$formula) <- 
      environment(object$formula)             
#    switch(offsetName, 
#           "1"     = { 
#                       if( (length(attr(Terms, "term.labels")) == 0) &&
#                           !is.null(attr(Terms, "offset")) )
#                       {
#                         ..offset <- modOff + offsetCoef[i]*.offset
#                          nf <- as.formula(".~. -1")
#                          nf <- update(.object$formula, nf)
##                         assign("..offset", ..offset, envir=sys.frame())
##                         assign("..offset", ..offset, pos=1)		## 20.05.13		
##                         nf <- as.formula(
##                                 paste(deparse(formula(object)[[2]],
##                                               width.cutoff = 500), 
##                                       " ~ -1 + offset(..offset)",
##                                       collapse = ""))
#                         .object$formula <- nf
#                         .object$offset <- ..offset
#                         environment(.object$formula) <- 
#                                  environment(object$formula)
#                       }
#                       else
#                       {
##                         nf <- as.formula(
##                                 paste(deparse(formula(object), 
##                                               width.cutoff = 500), 
##                                       "-1 + offset(", offsetCoef[i], 
##                                       "*.offset)", collapse = ""))
#                          nf <- as.formula(".~. -1")
#                          nf <- update(.object$formula, nf)
#                          .object$formula <- nf
#                          if(!is.null(.object$offset))
#                            .object$offset <- object$offset +   
#                                                offsetCoef[i]*.offset
#                          else
#                            .object$offset <- offsetCoef[i]*.offset  
#                          environment(.object$formula) <- 
#                                   environment(object$formula)
#                       }
#                     },
#           { 
#             nf <- as.formula(paste(".~.-", .offsetName))
#             nf <- update(.object$formula, nf)
##             nf <- as.formula(paste(deparse(nf, width.cutoff = 500), 
##                                    "+ offset(", offsetCoef[i], 
##                                    "* .offset)", collapse = ""))
#             .object$formula <- nf
#             if(!is.null(.object$offset))
#               .object$offset <- object$offset + offsetCoef[i]*.offset
#             else
#               .object$offset <- offsetCoef[i]*.offset  
#             environment(.object$formula) <- 
#                                  environment(object$formula)
#           })          
    .object <- eval(.object)
    summary.obj <- summary(.object)
    glmCDev[i] <- deviance(.object)
    if(!is.scalar)
#      glmCDet[i] <- det(summary.obj$cov.unscaled)
      glmCDet[i] <- determinant(summary.obj$cov.unscaled)$mod
  }
  if(trace) cat("\n")
  lp <- -1/2 * glmCDev
#  lmp <- if(!is.scalar) -1/2 * (glmCDev + log(glmCDet))
  lmp <- if(!is.scalar) -1/2 * (glmCDev + glmCDet)
         else lp
  lp <- spline(offsetCoef[is.finite(lp)], lp[is.finite(lp)], n)
  lmp <- spline(offsetCoef[is.finite(lmp)], lmp[is.finite(lmp)], n)
  cond.obj <- list(workspace = list(psi = offsetCoef,
                                    l.p = lp, l.mp = lmp))
  cond.obj$workspace$l.p$y <- ( cond.obj$workspace$l.p$y -
                                  max(cond.obj$workspace$l.p$y) )
  cond.obj$workspace$l.mp$y <- ( cond.obj$workspace$l.mp$y -
                                   max(cond.obj$workspace$l.mp$y) )
  if(!is.scalar)
  {
    s.mp <- predict(smooth.spline(lmp, all.knots=FALSE),
                    lmp$x, 1)
    glmCMLE <- predict(smooth.spline(s.mp$y, s.mp$x, all.knots=FALSE), 
                       0, 0)$y
    glmCSE <- sqrt(-predict(smooth.spline(s.mp$y, s.mp$x,
                                          all.knots=FALSE), 
                            0, 1)$y)
  }
  else
  {
    glmCMLE <- glmMLE
    glmCSE <- glmSE
  }
  Diff <- glmMLE - offsetCoef
  r.wald <- Diff/glmSE
  r.cwald <- (glmCMLE - offsetCoef)/glmCSE
  r.p <- sign(Diff) * sqrt(glmCDev - glmDeviance)
  Coef <- matrix(c(glmMLE, glmCMLE, glmSE, glmCSE), ncol=2)
  dimnames(Coef) <- list(c("uncond. ", "cond. "), 
                         c(" Estimate ", " Std. Error "))
##
##--> check for numerical instabilities
##
  condition <- abs(offsetCoef - glmMLE) > tms * glmSE
  poly.ord <- round(pts/2)
  aux.wmat1 <- matrix( rep(r.wald, poly.ord), ncol=poly.ord ) ^
                 matrix( rep(1:poly.ord, pts), ncol=poly.ord, 
                         byrow=TRUE )
  aux.wmat0 <- matrix( rep(r.wald, poly.ord), ncol=poly.ord ) ^
                 matrix( rep(0:(poly.ord-1), pts), ncol=poly.ord, 
                         byrow=TRUE)
  aux.pmat1 <- matrix( rep(r.p, poly.ord), ncol=poly.ord ) ^
                matrix( rep(1:poly.ord, pts), ncol=poly.ord, 
                        byrow=TRUE )
  aux.pmat0 <- matrix( rep(r.p, poly.ord), ncol=poly.ord ) ^
                matrix( rep(0:(poly.ord-1), pts), ncol=poly.ord, 
                        byrow=TRUE)
##
  rhoMax <- if(!is.scalar)
#               sqrt(glmDet * glmCDet) * glmSE
               sqrt(exp(glmDet + glmCDet)) * glmSE
             else rep(1, length=lengthOC)
  limINF <- predict(smooth.spline(lp, all.knots=FALSE),
                    glmMLE, 3)$y/6 * glmSE^3
  limNP <- if(!is.scalar)
#             predict(smooth.spline(offsetCoef, log(glmCDet), 
             predict(smooth.spline(offsetCoef, glmCDet, 
                                   all.knots=FALSE),
                     glmMLE, 1)$y/2 * glmSE
           else NULL
##
##-- L.-R.
##
  rhoLR <- rhoMax * r.wald
  rhoLRclow <- rhoLR * (exp(Diff) - 1)/Diff
  rhoLRcup <-  rhoLR * (1 - exp( - Diff))/Diff
  rhoLR <- ifelse(sign(r.p*rhoLR) < 0, -rhoLR, rhoLR)
  qT <- spline(x = offsetCoef[is.finite(rhoLR)], 
               y = rhoLR[is.finite(rhoLR)], n)  
  rhoLRclow <- ifelse(sign(r.p*rhoLRclow) < 0, -rhoLRclow, rhoLRclow)
  rhoLRcup <- ifelse(sign(r.p*rhoLRcup) < 0, -rhoLRcup, rhoLRcup)
##
  rhoLR.tmp <- (1/rhoLR - 1/r.p)
  aux.mod <- lm.fit(aux.pmat0[condition,], rhoLR.tmp[condition])
  rhoLR.tmp <- aux.pmat0%*%aux.mod$coef
  rhoLR.tmp <- ifelse(sign(r.p*rhoLR) < 0, -rhoLR.tmp, rhoLR.tmp)
  rhoLR <- ifelse(condition, 1/rhoLR - 1/r.p, rhoLR.tmp)
  LR <- pnorm(r.p) - dnorm(r.p) * rhoLR
##
  rhoLR.tmp <- (1/rhoLRclow - 1/r.p)
  aux.mod <- lm.fit(aux.pmat0[condition,], rhoLR.tmp[condition])
  rhoLR.tmp <- aux.pmat0%*%aux.mod$coef
  rhoLR.tmp <- ifelse(sign(r.p*rhoLRclow) < 0, -rhoLR.tmp, rhoLR.tmp)
  rhoLR <- ifelse(condition, 1/rhoLRclow - 1/r.p, rhoLR.tmp)
  LRclow <- pnorm(r.p) - dnorm(r.p) * rhoLR
##
  rhoLR.tmp <- (1/rhoLRcup - 1/r.p)
  aux.mod <- lm.fit(aux.pmat0[condition,], rhoLR.tmp[condition])
  rhoLR.tmp <- aux.pmat0%*%aux.mod$coef
  rhoLR.tmp <- ifelse(sign(r.p*rhoLRcup) < 0, -rhoLR.tmp, rhoLR.tmp)
  rhoLR <- ifelse( condition, 1/rhoLRcup - 1/r.p, rhoLR.tmp )
  LRcup <- pnorm(r.p) - dnorm(r.p) * rhoLR
##
##-- r*
##
  aux.mod <- lm.fit(aux.wmat1, r.p)
  logRho <- ifelse(condition, log(r.p/r.wald),
                              log(aux.wmat0%*%aux.mod$coef))
  rho <- logRho - log(abs(rhoMax))
##
  rho.tmp <- ifelse(sign(r.p*rho) < 0, -rho, rho)
  aux.mod <- lm.fit(aux.pmat1, rho.tmp)
  rho.tmp <- aux.pmat0%*%aux.mod$coef
  rho.tmp <- ifelse(sign(r.p*rho) < 0, -rho.tmp, rho.tmp)
  rho <- ifelse(condition, rho/r.p, rho.tmp)
##
  r.mp <- r.p - rho
  r.mp.clow <- r.mp + log((exp(Diff) - 1)/Diff)/r.p
  r.mp.cup <- r.mp + log((1 - exp( - Diff))/Diff)/r.p
##
##-- INF/NP
##
  inf <- inf.rp <- logRho
  np <- np.rp <- if(!is.scalar) -log(rhoMax)
                 else NULL
##
  inf.tmp <- ifelse(sign(r.p*inf) < 0, -inf, inf)
  aux.mod <- lm.fit(aux.pmat1, inf.tmp)
  inf.tmp <- aux.pmat0%*%aux.mod$coef
  inf.tmp <- ifelse(sign(r.p*inf) < 0, -inf.tmp, inf.tmp)
  inf <- inf.tmp
##
  if(!is.scalar)
  {
    np.tmp <- ifelse(sign(r.p*np) < 0, -np, np)
    aux.mod <- lm.fit(aux.pmat1, np.tmp)
    np.tmp <- aux.pmat0%*%aux.mod$coef
    np.tmp <- ifelse(sign(r.p*np) < 0, -np.tmp, np.tmp)
    np <- np.tmp
  }
##
##--> stop
##
  cond.obj$workspace <-
    c(cond.obj$workspace,
      list( r.e       = spline(offsetCoef[is.finite(r.wald)], 
                               r.wald[is.finite(r.wald)], n),
            r.ce      = spline(offsetCoef[is.finite(r.cwald)], 
                               r.cwald[is.finite(r.cwald)], n),
            r.p       = spline(offsetCoef[is.finite(r.p)], 
                               r.p[is.finite(r.p)], n), 
            r.mp      = spline(offsetCoef[is.finite(r.mp)], 
                               r.mp[is.finite(r.mp)], n),
            r.mp.clow = spline(offsetCoef[is.finite(r.mp.clow)], 
                               r.mp.clow[is.finite(r.mp.clow)], n),
            r.mp.cup  = spline(offsetCoef[is.finite(r.mp.cup)], 
                               r.mp.cup[is.finite(r.mp.cup)], n) ))
  lr <- spline(offsetCoef[is.finite(LR)], LR[is.finite(LR)], n)
  lr$y[lr$y < 0] <- 0
  lr$y[lr$y > 1] <- 1
  lr.clow <- spline(offsetCoef[is.finite(LRclow)], 
                    LRclow[is.finite(LRclow)], n)
  lr.clow$y[lr.clow$y < 0] <- 0
  lr.clow$y[lr.clow$y > 1] <- 1
  lr.cup <-  spline(offsetCoef[is.finite(LRcup)], 
                    LRcup[is.finite(LRcup)], n)
  lr.cup$y[lr.cup$y < 0] <- 0
  lr.cup$y[lr.cup$y > 1] <- 1
  cond.obj$workspace <-
    c(cond.obj$workspace,
      list( lr        = lr, 
            lr.clow   = lr.clow,
            lr.cup    = lr.cup,
            q         = qT, 
            inf       = spline(offsetCoef[is.finite(inf)], 
                               inf[is.finite(inf)], n),
            inf.rp    = spline(offsetCoef[is.finite(inf.rp)], 
                               inf.rp[is.finite(inf.rp)], n),
            limINF    = list(x = glmMLE, y = limINF),
            np        = if(!is.scalar) 
                          spline(offsetCoef[is.finite(np)], 
                                 np[is.finite(np)], n)
                        else NULL,
            np.rp     = if(!is.scalar) 
                          spline(offsetCoef[is.finite(np.rp)], 
                                 np.rp[is.finite(np.rp)], n)
                        else NULL,
            limNP     = list(x = glmMLE, y = limNP) ))
  cond.obj <- c(cond.obj, 
                list( coefficients = Coef,
                      call         = m,
                      formula      = formula,
                      family       = as.name(object$family[[1]]),
                      offset       = as.name(offsetName),
                      diagnostics  = c(INF = 
                               max(abs(cond.obj$workspace$inf$y)),
                                       NP  = if(!is.scalar)
                               max(abs(cond.obj$workspace$np$y))
                                             else NULL),
                      n.approx     = pts,
                      omitted.val  = c(glmMLE - tms * glmSE,
                                       glmMLE + tms * glmSE),
                      is.scalar    = is.scalar ))
  attr(cond.obj, "class") <- "cond"
  cond.obj
}

summary.cond <- function(object, alpha = 0.05, test = NULL, 
                         all = FALSE, coef = TRUE, 
                         int = ifelse((is.null(test) || all), 
                                 TRUE, FALSE), 
                         digits = NULL, ...)
{
  m <- match.call()
  dim.alpha <- length(alpha)
  if( !is.null(test) )
    dim.test <- length(test)
  alpha.quant <- c(qnorm(1 - alpha/2), qnorm(1 - alpha))
#  attach(object$workspace, warn.conflicts = FALSE)
#  on.exit( detach() )
  r.p <- object$workspace$r.p
  r.mp <- object$workspace$r.mp
  r.mp.cup <- object$workspace$r.mp.cup
  r.mp.clow <- object$workspace$r.mp.clow
  lr <- object$workspace$lr
  lr.cup <- object$workspace$lr.cup
  lr.clow <- object$workspace$lr.clow
  cf <- object$coefficients
  is.scalar <- object$is.scalar
  lr.ss <- try(smooth.spline(lr$y, lr$x), silent=TRUE)
  lr.try <- !inherits(lr.ss, "try-error")
  if(int)
  {
    CI.mle <- cf[1,1] - cf[1,2] * c(alpha.quant, -alpha.quant)
    CI.cmle <- cf[2,1] - cf[2,2] * c(alpha.quant, -alpha.quant)
    CI.rp <- predict(smooth.spline(r.p$y, r.p$x),
                     c(alpha.quant, -alpha.quant), 0)$y
    CI.rmp <- predict(smooth.spline(r.mp$y, r.mp$x),
                      c(alpha.quant, -alpha.quant), 0)$y
    CI.rmp.corr <- c(predict(smooth.spline(r.mp.cup$y, r.mp.cup$x),
                             alpha.quant, 0)$y,
                     predict(smooth.spline(r.mp.clow$y, r.mp.clow$x),
                             -alpha.quant, 0)$y)
    if( lr.try)
    {
      CI.lr <- predict(lr.ss,
                       c(1-alpha/2, 1-alpha, alpha/2, alpha), 0)$y
      CI.lr.corr <- c(predict(smooth.spline(lr.cup$y, lr.cup$x),
                               c(1-alpha/2, 1-alpha), 0)$y,
                      predict(smooth.spline(lr.clow$y, lr.clow$x),
                              c(1-alpha/2, 1-alpha), 0)$y)
    }
    else
      CI.lr <- CI.lr.corr <- rep(NA, length(CI.mle))
    CI <- t(matrix(c(CI.mle, CI.cmle, CI.rp, CI.rmp, CI.rmp.corr,
                     CI.lr, CI.lr.corr),
                   ncol = 4*dim.alpha, byrow = TRUE))
    dimnames(CI) <- list(c(paste("Lower conf. bound (2sd,", 
                                 100*(1-alpha), "%)"),
                           paste("Lower conf. bound (1sd,", 
                                 100*(1-alpha), "%)"),
                           paste("Upper conf. bound (2sd,", 
                                 100*(1-alpha), "%)"),
                           paste("Upper conf. bound (1sd,", 
                                 100*(1-alpha), "%)")),
                 c("Wald pivot                             ",
                   "Wald pivot (cond. MLE)                 ",
                   "Likelihood root                        ",
                   "Modified likelihood root               ",
                   "Modified likelihood root (cont. corr.) ",
                   "Lugannani-Rice approx.                 ",
                   "Lugannani-Rice approx. (cont. corr.)   "))
  }
  old.test <- test
  if(!is.null(test))
  {
    test.mle <- (cf[1,1] - test) / cf[1,2]
    test.cmle <- (cf[2,1] - test) / cf[2,2]
    test.rp <- predict(smooth.spline(r.p), test, 0)$y
    test.rmp <- predict(smooth.spline(r.mp), test, 0)$y
    test.rmp.corr <- ifelse(test.rmp < 0,
                          predict(smooth.spline(r.mp.clow), test, 0)$y,
                          predict(smooth.spline(r.mp.cup), test, 0)$y)
    if( lr.try )
    {
      test.lr <- predict(smooth.spline(lr), test, 0)$y
      test.lr.corr <- ifelse(test.lr < 0.5,
                              predict(smooth.spline(lr.clow), test, 0)$y,
                            1 - predict(smooth.spline(lr.cup), test, 0)$y)
      test.lr <- ifelse(test.lr > 0.5, 1 - test.lr, test.lr)
    }
    else
      test.lr <- test.lr.corr <- rep(NA, dim.test)
    TEST <- t(matrix(c(test.mle,
                       ifelse(test.mle > 0, (1 - pnorm(test.mle)),
                                            pnorm(test.mle)),
                       test.cmle,
                       ifelse(test.cmle > 0, (1-pnorm(test.cmle)),
                                             pnorm(test.cmle)),
                       test.rp,
                       ifelse(test.rp > 0, (1-pnorm(test.rp)),
                                           pnorm(test.rp)),
                       test.rmp,
                       ifelse(test.rmp > 0, (1-pnorm(test.rmp)),
                                            pnorm(test.rmp)),
                       test.rmp.corr,
                     ifelse(test.rmp.corr > 0, (1-pnorm(test.rmp.corr)),
                                               pnorm(test.rmp.corr)),
                       rep(NA, dim.test), test.lr,
                       rep(NA, dim.test), test.lr.corr),
                    nrow = 2 * dim.test))
    dimnames(TEST) <- list(
                 c("Wald pivot                             ",
                   "Wald pivot (cond. MLE)                 ",
                   "Likelihood root                        ",
                   "Modified likelihood root               ",
                   "Modified likelihood root (cont. corr.) ",
                   "Lugannani-Rice approx.                 ",
                   "Lugannani-Rice approx. (cont. corr.)   "),
                 c(paste("statistic (", object$offset, 
                          "=", test, ")  "),
                   paste("tail prob. (", object$offset, 
                          "=", test, ")        ")))
   qT <- object$workspace$q
   qT <- predict(smooth.spline(qT, all.knots=FALSE), test, 0)$y
   qT <- matrix(qT, ncol=1)
   dimnames(qT) <- list(paste("(", object$offset, "=", test, ")"),
                        "q")
   TEST <- list(stats=TEST, qTerm=qT)
  }
  summary.object <- list( coefficients = cf,
                          conf.int = if(int) CI,
                          signif.tests = if(!is.null(test)) TEST ,
                          call = object$call,
                          formula = object$formula,
                          family = object$family,
                          offset = object$offset,
                          alpha = alpha,
                          hypotheses = test,
                          diagnostics = object$diagnostics,
                          n.approx = object$n.approx,
                          all = all, cf = coef, int = int,
                          is.scalar = is.scalar, digits = digits )
  class(summary.object) <- "summary.cond"
  summary.object
}

plot.cond <- function(x = stop("nothing to plot"), from = x.axis[1], 
                      to = x.axis[n], which = NULL, alpha = 0.05, 
                      add.leg = TRUE, loc.leg = FALSE, add.labs = TRUE, 
                      cex = 0.7, cex.lab = 1, cex.axis = 1, cex.main = 1, 
                      lwd1 = 1, lwd2 = 2, lty1 = "solid", lty2 = "dashed", 
                      col1 = "black", col2 = "blue", tck = 0.02, las = 1, 
                      adj = 0.5, lab = c(15, 15, 5), ...)
{
  offset <- x$offset
  if( offset == "1" )
    offset <- "Intercept"
  is.scalar <- x$is.scalar
  choices <- c("All",
               ifelse(!is.scalar,
                 "Profile and modified profile log likelihoods",
                 "Profile log likelihood"),
               ifelse(!is.scalar,
                 "Profile and modified profile likelihood ratios",
                 "Profile likelihood ratio"),
               "Profile and modified likelihood roots",
               "Modified and continuity corrected likelihood roots",
               "Lugannani-Rice approximations",
               "Confidence intervals",
               ifelse(!is.scalar,
                 "Diagnostics based on INF/NP decomposition\n",
                 "Diagnostics based on INF decomposition\n"))
  tmenu <- paste("", choices)
  if( is.null(which) )
    pick <- menu(tmenu, 
                 title="\n Make a plot selection (or 0 to exit)\n")
  else if( !match(which, 2:8, nomatch=FALSE) )
          stop("choice not valid")
       else pick <- which
  if(pick == 0)
    stop(" no graph required ! ")
#  attach(x$workspace, warn.conflicts = FALSE)
#  on.exit(invisible(detach()))
  r.p <- x$workspace$r.p
  r.mp <- x$workspace$r.mp
  r.mp.cup <- x$workspace$r.mp.cup
  r.mp.clow <- x$workspace$r.mp.clow
  l.p <- x$workspace$l.p
  l.mp <- x$workspace$l.mp
  r.e <- x$workspace$r.e
  r.ce <- x$workspace$r.ce
  lr.clow <- x$workspace$lr.clow
  lr.cup <- x$workspace$lr.cup
  lr <- x$workspace$lr
  inf <- x$workspace$inf
  limINF <- x$workspace$limINF
  inf.rp <- x$workspace$inf.rp
  np <- x$workspace$np
  limNP <- x$workspace$limNP
  np.rp <- x$workspace$np.rp
  coeff <- x$coefficients
  x.alpha <- qnorm(1 - alpha/2)
  CI.rp <- predict(smooth.spline(r.p$y, r.p$x), 
                   c(x.alpha, -x.alpha), 0)$y
  CI.rmp <- predict(smooth.spline(r.mp$y, r.mp$x), 
                    c(x.alpha, -x.alpha), 0)$y
  CI.rmp.corr <- c(predict(smooth.spline(r.mp.cup$y,r.mp.cup$x), 
                           x.alpha, 0)$y,
                   predict(smooth.spline(r.mp.clow$y, r.mp.clow$x), 
                           -x.alpha, 0)$y)
  CI.mle <- (coeff[1,1] - coeff[1,2]*c(x.alpha, -x.alpha))
  CI.cmle <- (coeff[2,1] - coeff[2,2]*c(x.alpha, -x.alpha))
  from.ic <- min(c(CI.rp[1], CI.rmp[1], CI.rmp.corr[1], CI.mle[1], 
                   CI.cmle[1]))
  to.ic <- max(c(CI.rp[2], CI.rmp[2], CI.rmp.corr[2], CI.mle[2], 
               CI.cmle[2]))
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par), add=TRUE)
  repeat
  {
    invisible(par(tck=tck, las=las, adj=adj, lab=lab, ...))
    switch(pick,
		"1" = {
			invisible(split.screen(figs = c(1,2)))
##               ---------------
           screen(1)
           n <- length(l.p$x)
           x.axis <- l.p$x
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ((l.p$y > -3) | (l.mp$y > -3))
           plot(x.axis[condition], l.p$y[condition], 
                type = "n", lty = lty1, lwd = lwd1, 
                ylim=c(max(-3, 
                       min(l.p$y[condition], l.mp$y[condition])), 0),
                xlab = "", ylab = "", cex.lab = cex.lab, 
                cex.axis = cex.axis, cex.main = cex.main, ...)
           lines(x.axis[condition], l.p$y[condition], lty = lty1, 
                 col = col1, lwd=lwd1, ...)
           if( !is.scalar )
             lines(x.axis[condition], l.mp$y[condition], lty=lty1, 
                   col=col2, lwd=lwd2, ...)
           if(add.labs)
             title(xlab = paste("coefficient of", offset),
                   ylab = "log likelihood",
                   main = ifelse(!is.scalar,
                                 "Profile and modified profile log likelihoods",
                                 "Profile log likelihood"),
                   ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- (x.axis[l.mp$y == max(l.mp$y)] +
                         x.axis[condition][1])/2
               y.lim <- max(-3,
                            min(l.p$y[condition], l.mp$y[condition]))/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile log likelihood",
                        "modified profile log likelihood"), cex=cex,
                      lty=lty1, lwd=c(lwd1,lwd2), bty="n", 
                      col=c(col1,col2), ...)
             else 
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      "profile log likelihood", cex=cex,
                      lty=lty1, lwd=lwd1, bty="n", col=col1,
                      ...)
           }
##               ---------------
           screen(2)
           conf.limit <-  qchisq(1-alpha, 1)
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.e$y)^2 < (conf.limit + 1)) |
                          ((r.ce$y)^2 < (conf.limit + 1)) |
                          (-2 * l.p$y < (conf.limit + 1)) |
                          (-2 * l.mp$y < (conf.limit + 1)) )
           if( !is.scalar )
           {
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            ((r.ce$y)^2)[condition],
                            (-2 * l.p$y)[condition] ),
                     xlab="", ylab="", type="l", lty=c(lty2,lty2,lty1), 
                     lwd=c(lwd1,lwd2,lwd1), col=c(col1,col2,col1), 
                     cex.lab=cex.lab, cex.main=cex.main, 
                     cex.axis=cex.axis, ...)
             lines(x.axis[condition], (-2 * l.mp$y)[condition], 
                   lty=lty1, lwd=lwd2, col=col2, ...)
           }
           else
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            (-2 * l.p$y)[condition] ),
                     xlab="", ylab="", type="l", lty=c(lty2,lty1), 
                     lwd=c(lwd1,lwd1), col=c(col1,col1), 
                     cex.lab=cex.lab, cex.main=cex.main, 
                     cex.axis=cex.axis, ...)
           abline(h=conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab = paste("coefficient of", offset),
                   ylab = "likelihood ratio",
                   main = ifelse(!is.scalar,
                                 "Profile and modified profile likelihodo ratios",
                                 "Profile likelihood ratio"),
                   ...)
 	   if(add.leg)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- ( x.axis[l.p$y == max(l.p$y)] +
                          x.axis[condition][1] )/2
               tt <- max(c(max(((r.e$y)^2)[condition]),
                           max(((r.ce$y)^2)[condition]),
                           max((-2 * l.p$y)[condition]),
                           max((-2 * l.mp$y)[condition])))
               y.lim <- if(tt > conf.limit) (tt+conf.limit)/2
                        else tt/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile likelihood ratio", "modified profile LR",
                        "Wald pivot",
                        "Wald pivot (cond. MLE)"),
                      lty=c(lty1,lty1,lty2,lty2), 
                      lwd=c(lwd1,lwd2,lwd1,lwd2),
                      col=c(col1,col2,col1,col2), bty="n", 
                      cex=cex, ...)
             else
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile LR", "Wald pivot"),
                      lty=c(lty1,lty2), lwd=c(lwd1,lwd1),
                      col=c(col1,col1), bty="n", cex=cex, ...)

           }
##               ---------------
                        invisible(close.screen(all.screens=TRUE))
                        par(ask=TRUE)
                        invisible(split.screen(figs=c(2,2)))
##               ---------------
           conf.limit <- qnorm(1 - alpha/2)
           screen(1)
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.e$y < 4) & (r.e$y > -4)) |
                          ((r.ce$y < 4) & (r.ce$y > -4)) |
                          ((r.p$y < 4) & (r.p$y > -4)) |
                          ((r.mp$y < 4) & (r.mp$y > -4)) )
           plot(0, 0, type="n", xlim=range(x.axis[condition]), ylim=c(-4, 4),
                xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.e$y < 4) & (r.e$y > -4)
	   lines(x.axis[condition], r.e$y[condition], lty=lty2, col=col1,
                 lwd=lwd1, ...)
           if( !is.scalar )
           {
#             condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.ce$y < 4) & (r.ce$y > -4)
             lines(x.axis[condition], r.ce$y[condition], lty=lty2, 
                   col=col2, lwd=lwd2, ...)
           }
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.p$y < 4) & (r.p$y > -4)
	   lines(x.axis[condition], r.p$y[condition], lty=lty1, 
                 lwd=lwd1, col=col1, ...)
           n <- length(r.mp$x)
           x.axis <- r.mp$x
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp$y < 4) & (r.mp$y > -4)
     	   lines(x.axis[condition], r.mp$y[condition], lty=lty1, 
                 lwd=lwd2, col=col2, ... )
    	   abline(h=0, lty="dotted", ...)
	   abline(h=conf.limit, lty="dotted", ...)
	   abline(h=-conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab=paste("coefficient of", offset),
                   ylab="likelihood root",
                   main="Profile and modified likelihood roots",
                   ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- x.axis[condition][2]
               y.lim <- -2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root",
                        "modified likelihood root",
                        "Wald pivot",
                        "Wald pivot (cond. MLE)"),
                      lty=c(lty1,lty1,lty2, lty2), 
                      lwd=c(lwd1,lwd2,lwd1,lwd2), 
                      bty="n", cex=cex, 
                      col=c(col1, col2, col1, col2), ...)
             else
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root",
                        "modified likelihood root",
                        "Wald pivot"),
                      lty=c(lty1,lty1,lty2), 
                      lwd=c(lwd1,lwd2,lwd1), 
                      bty="n", cex=cex, 
                      col=c(col1, col2, col1), ...)
           }
##               ---------------
           screen(2)
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.mp$y < 4) & (r.mp$y > -4)) |
                          ((r.mp.clow$y < 4) & (r.mp.clow$y > -4)) |
                          ((r.mp.cup$y < 4) & (r.mp.cup$y > -4)) )
           plot(0, 0, type="n", xlim=range(x.axis[condition]), ylim=c(-4, 4),
                xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp$y < 4) & (r.mp$y > -4)
           lines(x.axis[condition], r.mp$y[condition], lty=lty1, 
                 lwd=lwd1, col=col1, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp.clow$y < 4) & (r.mp.clow$y > -4)
           lines(x.axis[condition], r.mp.clow$y[condition], lty=lty1, 
                 lwd=lwd1, col=col2, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp.cup$y < 4) & (r.mp.cup$y > -4)
           lines(x.axis[condition], r.mp.cup$y[condition], lty=lty1, 
                 lwd=lwd1, col=col2, ...)
	   abline(h=0, lty="dotted", ...)
	   abline(h=conf.limit, lty="dotted", ...)
	   abline(h=-conf.limit, lty="dotted", ...)
           if(add.labs)
              title(xlab = paste("coefficient of", offset),
                    ylab = "likelihood root",
                    main = "Modified likelihood roots",
                    ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- x.axis[condition][2]
               y.lim <- -2
             }
	     legend(c(x.lim, x.lim), c(y.lim, y.lim),
                    c("modified likelihood root",
	   	      "modified likelihood root \n (cont. corr.)"),
                    lty=c(lty1,lty1) , lwd=c(lwd1,lwd2), bty="n",
                    cex=cex, col=c(col1,col2), ...)
           }
##               ---------------
           screen(3)
           n <- length(r.mp$x)
           x.axis <- r.mp$x
           condition <- (x.axis >= from) & (x.axis <= to)
           matplot(x.axis[condition],
                   cbind(lr.clow$y[condition], lr.cup$y[condition]),
                   xlab="", ylab= "", type="l", lty=lty1, lwd=lwd1,
                   ylim=c(0,1), col=col2, cex.lab=cex.lab, 
                   cex.axis=cex.axis, cex.main=cex.main, ...)
           lines(x.axis[condition], lr$y[condition], lty=lty1, lwd=lwd1, 
                 col=col1, ...)
           if(add.labs)
              title(xlab = paste("coefficient of", offset),
                    ylab = "tail area approximation",
                    main = "Lugannani-Rice approximations",
                    ...)
	   if(add.leg == TRUE)
           {
                   if(loc.leg)
                   {
                       cat("Choose legend position\n")
                       loci <- locator(1)
                       x.lim <- loci$x
                       y.lim <- loci$y
                   }
                   else
                   {
                       tt <- lr$y[condition][1]
                       x.lim <- if(tt < 0.5)
                                  x.axis[condition][
                                      (trunc(length(x.axis[condition])/2))]
                                else x.axis[condition][2]
                       y.lim <- if(tt < 0.5) (tt+1)/2 else tt/4
                   }
	   	legend(c(x.lim, x.lim), c(y.lim, y.lim),
                          c("tail area approximation",
                            "tail area approx. \n (cont. corr.)"),
                          lty=c(lty1,lty1), lwd=c(lwd1,lwd2), bty="n",
                          cex=cex, col=c(col1,col2), ...)
           }
##               ---------------
           screen(4)
           invisible(par(lab=c(30, 5, 5)))
           plot(coeff[1,1], 2.5, type="n", xlim=c(from.ic, to.ic), 
                ylim=c(0, ifelse(!is.scalar, 6, 5)), xlab="", 
                ylab="", cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, ...)
           segments(CI.rmp.corr[1], 1, CI.rmp.corr[2], 1, lty=lty1,
                    lwd=lwd2, col=col2, ...)
           segments(CI.rmp[1], 2, CI.rmp[2], 2, lty=lty1, lwd=lwd2, 
                    col=col2, ...)
           segments(CI.rp[1], 3, CI.rp[2], 3, lty=lty1, lwd=lwd2, 
                    col=col1, ...)
           if( !is.scalar )
           {
             segments(CI.cmle[1], 4, CI.cmle[2], 4, lty=lty1, lwd=lwd2, 
                      col=col2, ...)
             segments(CI.mle[1], 5, CI.mle[2], 5, lty=lty1, lwd=lwd2, 
                      col=col1, ...)
           }
           else
             segments(CI.mle[1], 4, CI.mle[2], 4, lty=lty1, lwd=lwd2, 
                      col=col1, ...)
	   if(add.labs)
              title(xlab=paste("coefficient of", offset),
                    ylab="confidence interval",
                    main=paste(100*(1-alpha), 
                               "% Confidence Intervals "),
                    ...)
           if(add.leg == TRUE)
           {
             text(coeff[2,1], 1.2,
                  "modified likelihood root (cont.corr.)", cex=cex, ...)
             text(coeff[2,1], 2.2, "modified likelihood root",
                  cex=cex, ...)
             text(coeff[2,1], 3.2, "likelihood root", cex=cex, ...)
             if( !is.scalar )
             {
               text(coeff[2,1], 4.2, "Wald pivot (cond. MLE)", 
                    cex=cex, ...)
               text(coeff[2,1], 5.2, "Wald pivot",
                     cex=cex, ...)
             }
             else
               text(coeff[2,1], 4.2, "Wald pivot",
                     cex=cex, ...)               
           }
##               ---------------
			invisible(close.screen(all.screens=TRUE))
                        par(ask=TRUE)
                        par(pty="s")
                        if(is.scalar) split.screen(figs=c(1,2))
                        else split.screen(figs=c(2,2))
##               ---------------
           screen(1)
           plot(inf, type="l", xlab="", ylab="", lty=lty1, 
                lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, col=col1, ...)
           points(limINF, pch="o", cex=cex.axis, col=col1, ...)
           abline(h=0.2, lty="dotted", ...)
           abline(h=-0.2, lty="dotted", ...)
           if(add.labs)
             title(xlab=paste("coefficient of", offset), ylab="INF",
                   main="INF correction term", ...)
           screen(2)
           plot(inf.rp, type="l", xlab="", ylab="", lty=lty1, 
                   lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, col=col1, ...)
           if(add.labs)
             title(xlab="likelihood root", ylab="INF",
                   main="INF correction term", ...)
           if( !is.scalar )
           {
             screen(3)
             plot(np, type ="l", xlab="", ylab="", lty=lty1, 
                  lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                  cex.main=cex.main, col=col1, ...)
             points(limNP, pch="o", cex=cex.axis, col=col1, ...)
             abline(h=0.2, lty=3, ...)
             abline(h=-0.2, lty=3, ...)
             if(add.labs)
               title(xlab=paste("coefficient of", offset), ylab="NP",
                     main="NP correction term", ...)
             screen(4)
             plot(np.rp, type="l", xlab="", ylab="", lty=lty1, 
                  lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                  cex.main=cex.main, col=col1, ...)
             if(add.labs)
               title(xlab="likelihood root", ylab="NP",
                     main="NP correction term", ...)
           }
##               ---------------
                        close.screen(all.screens=TRUE)
                        par(pty="m")
                        par(ask=FALSE)
		}
		,
		"2"={
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           n <- length(l.p$x)
           x.axis <- l.p$x
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ((l.p$y > -3) | (l.mp$y > -3))
           plot(x.axis[condition], l.p$y[condition], 
                type = "n", lty = lty1, lwd = lwd1, 
                ylim=c(max(-3, 
                       min(l.p$y[condition], l.mp$y[condition])), 0),
                xlab = "", ylab = "", cex.lab = cex.lab, 
                cex.axis = cex.axis, cex.main = cex.main, ...)
           lines(x.axis[condition], l.p$y[condition], lty = lty1, 
                 col = col1, lwd=lwd1, ...)
           if( !is.scalar )
             lines(x.axis[condition], l.mp$y[condition], lty=lty1, 
                   col=col2, lwd=lwd2, ...)
           if(add.labs)
             title(xlab = paste("coefficient of", offset),
                   ylab = "log likelihood",
                   main = ifelse(!is.scalar,
                                 "Profile and modified profile log likelihoods",
                                 "Profile log likelihood"),
                   ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- (x.axis[l.mp$y == max(l.mp$y)] +
                         x.axis[condition][1])/2
               y.lim <- max(-3,
                            min(l.p$y[condition], l.mp$y[condition]))/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile log likelihood",
                        "modified profile log likelihood"), cex=cex,
                      lty=lty1, lwd=c(lwd1,lwd2), bty="n", 
                      col=c(col1,col2), ...)
             else 
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      "profile log likelihood", cex=cex,
                      lty=lty1, lwd=lwd1, bty="n", col=col1,
                      ...)
           }
##               ---------------
                }
		,
		"3"={
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           n <- length(l.p$x)
           x.axis <- l.p$x
           conf.limit <-  qchisq(1-alpha, 1)
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.e$y)^2 < (conf.limit + 1)) |
                          ((r.ce$y)^2 < (conf.limit + 1)) |
                          (-2 * l.p$y < (conf.limit + 1)) |
                          (-2 * l.mp$y < (conf.limit + 1)) )
           if( !is.scalar )
           {
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            ((r.ce$y)^2)[condition],
                            (-2 * l.p$y)[condition] ),
                     xlab="", ylab="", type="l", lty=c(lty2,lty2,lty1), 
                     lwd=c(lwd1,lwd2,lwd1), col=c(col1,col2,col1), 
                     cex.lab=cex.lab, cex.main=cex.main, 
                     cex.axis=cex.axis, ...)
             lines(x.axis[condition], (-2 * l.mp$y)[condition], 
                   lty=lty1, lwd=lwd2, col=col2, ...)
           }
           else
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            (-2 * l.p$y)[condition] ),
                     xlab="", ylab="", type="l", lty=c(lty2,lty1), 
                     lwd=c(lwd1,lwd1), col=c(col1,col1), 
                     cex.lab=cex.lab, cex.main=cex.main, 
                     cex.axis=cex.axis, ...)
           abline(h=conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab = paste("coefficient of", offset),
                   ylab = "likelihood ratio",
                   main = ifelse(!is.scalar,
                                 "Profile and modified profile likelihood ratios",
                                 "Profile likelihood ratio"),
                   ...)
 	   if(add.leg)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- ( x.axis[l.p$y == max(l.p$y)] +
                          x.axis[condition][1] )/2
               tt <- max(c(max(((r.e$y)^2)[condition]),
                           max(((r.ce$y)^2)[condition]),
                           max((-2 * l.p$y)[condition]),
                           max((-2 * l.mp$y)[condition])))
               y.lim <- if(tt > conf.limit) (tt+conf.limit)/2
                        else tt/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile likelihood ratio", "modified profile LR",
                        "Wald pivot",
                        "Wald pivot (cond. MLE)"),
                      lty=c(lty1,lty1,lty2,lty2), 
                      lwd=c(lwd1,lwd2,lwd1,lwd2),
                      col=c(col1,col2,col1,col2), bty="n", 
                      cex=cex, ...)
             else
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile LR", "Wald pivot"),
                      lty=c(lty1,lty2), lwd=c(lwd1,lwd1),
                      col=c(col1,col1), bty="n", cex=cex, ...)

           }
##               ---------------
       	        }
		,
		"4"={
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           conf.limit <- qnorm(1 - alpha/2)
           n <- length(l.p$x)
           x.axis <- l.p$x
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.e$y < 4) & (r.e$y > -4)) |
                          ((r.ce$y < 4) & (r.ce$y > -4)) |
                          ((r.p$y < 4) & (r.p$y > -4)) |
                          ((r.mp$y < 4) & (r.mp$y > -4)) )
           plot(0, 0, type="n", xlim=range(x.axis[condition]), ylim=c(-4, 4),
                xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.e$y < 4) & (r.e$y > -4)
	   lines(x.axis[condition], r.e$y[condition], lty=lty2, col=col1,
                 lwd=lwd1, ...)
           if( !is.scalar )
           {
#             condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.ce$y < 4) & (r.ce$y > -4)
             lines(x.axis[condition], r.ce$y[condition], lty=lty2, 
                   col=col2, lwd=lwd2, ...)
           }
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.p$y < 4) & (r.p$y > -4)
	   lines(x.axis[condition], r.p$y[condition], lty=lty1, 
                 lwd=lwd1, col=col1, ...)
           n <- length(r.mp$x)
           x.axis <- r.mp$x
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp$y < 4) & (r.mp$y > -4)
     	   lines(x.axis[condition], r.mp$y[condition], lty=lty1, 
                 lwd=lwd2, col=col2, ... )
    	   abline(h=0, lty="dotted", ...)
	   abline(h=conf.limit, lty="dotted", ...)
	   abline(h=-conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab=paste("coefficient of", offset),
                   ylab="likelihood root",
                   main="Profile and modified likelihood root",
                   ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- x.axis[condition][2]
               y.lim <- -2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root",
                        "modified likelihood root",
                        "Wald pivot",
                        "Wald pivot (cond. MLE)"),
                      lty=c(lty1,lty1,lty2, lty2), 
                      lwd=c(lwd1,lwd2,lwd1,lwd2), 
                      bty="n", cex=cex, 
                      col=c(col1, col2, col1, col2), ...)
             else
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root",
                        "modified likelihood root",
                        "Wald pivot"),
                      lty=c(lty1,lty1,lty2), 
                      lwd=c(lwd1,lwd2,lwd1), 
                      bty="n", cex=cex, 
                      col=c(col1, col2, col1), ...)
           }
##               ---------------
		}
		,
		"5"={
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           conf.limit <- qnorm(1 - alpha/2)
           n <- length(r.mp$x)
           x.axis <- r.mp$x
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.mp$y < 4) & (r.mp$y > -4)) |
                          ((r.mp.clow$y < 4) & (r.mp.clow$y > -4)) |
                          ((r.mp.cup$y < 4) & (r.mp.cup$y > -4)) )
          plot(0, 0, type="n", xlim=range(x.axis[condition]), ylim=c(-4, 4),
                xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp$y < 4) & (r.mp$y > -4)
           lines(x.axis[condition], r.mp$y[condition], lty=lty1, 
                 lwd=lwd1, col=col1, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp.clow$y < 4) & (r.mp.clow$y > -4)
           lines(x.axis[condition], r.mp.clow$y[condition], lty=lty1, 
                 lwd=lwd1, col=col2, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                        (r.mp.cup$y < 4) & (r.mp.cup$y > -4)
           lines(x.axis[condition], r.mp.cup$y[condition], lty=lty1, 
                 lwd=lwd1, col=col2, ...)
	   abline(h=0, lty="dotted", ...)
	   abline(h=conf.limit, lty="dotted", ...)
	   abline(h=-conf.limit, lty="dotted", ...)
           if(add.labs)
              title(xlab = paste("coefficient of", offset),
                    ylab = "likelihood root",
                    main = "Modified likelihood root",
                    ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- x.axis[condition][2]
               y.lim <- -2
             }
	     legend(c(x.lim, x.lim), c(y.lim, y.lim),
                    c("modified likelihood root",
	   	      "modified likelihood root \n (cont. corr.)"),
                    lty=c(lty1,lty1) , lwd=c(lwd1,lwd2), bty="n",
                    cex=cex, col=c(col1,col2), ...)
           }
##               ---------------
		}
                ,
                "6"={
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           n <- length(r.mp$x)
           x.axis <- r.mp$x
           condition <- (x.axis >= from) & (x.axis <= to)
           matplot(x.axis[condition],
                   cbind(lr.clow$y[condition], lr.cup$y[condition]),
                   xlab="", ylab= "", type="l", lty=lty1, lwd=lwd1,
                   ylim=c(0,1), col=col2, cex.lab=cex.lab, 
                   cex.axis=cex.axis, cex.main=cex.main, ...)
           lines(x.axis[condition], lr$y[condition], lty=lty1, lwd=lwd1, 
                 col=col1, ...)
           if(add.labs)
              title(xlab = paste("coefficient of", offset),
                    ylab = "tail area approximation",
                    main = "Lugannani-Rice approximations",
                    ...)
	   if(add.leg == TRUE)
           {
                   if(loc.leg)
                   {
                       cat("Choose legend position\n")
                       loci <- locator(1)
                       x.lim <- loci$x
                       y.lim <- loci$y
                   }
                   else
                   {
                       tt <- lr$y[condition][1]
                       x.lim <- if(tt < 0.5)
                                  x.axis[condition][
                                      (trunc(length(x.axis[condition])/2))]
                                else x.axis[condition][2]
                       y.lim <- if(tt < 0.5) (tt+1)/2 else tt/4
                   }
	   	legend(c(x.lim, x.lim), c(y.lim, y.lim),
                          c("tail area approximation",
                            "tail area approx. \n (cont. corr.)"),
                          lty=c(lty1,lty1), lwd=c(lwd1,lwd2), bty="n",
                          cex=cex, col=c(col1,col2), ...)
           }
##               ---------------
                },
                "7"={
                        invisible(par(lab=c(30, 5, 5), mfrow=c(1,1)))
##               ---------------
           conf.limit <- qnorm(1 - alpha/2)
           plot(coeff[1,1], 2.5, type="n", xlim=c(from.ic, to.ic), 
                ylim=c(0, ifelse(!is.scalar, 6, 5)), xlab="", 
                ylab="", cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, ...)
           segments(CI.rmp.corr[1], 1, CI.rmp.corr[2], 1, lty=lty1,
                    lwd=lwd2, col=col2, ...)
           segments(CI.rmp[1], 2, CI.rmp[2], 2, lty=lty1, lwd=lwd2, 
                    col=col2, ...)
           segments(CI.rp[1], 3, CI.rp[2], 3, lty=lty1, lwd=lwd2, 
                    col=col1, ...)
           if( !is.scalar )
           {
             segments(CI.cmle[1], 4, CI.cmle[2], 4, lty=lty1, lwd=lwd2, 
                      col=col2, ...)
             segments(CI.mle[1], 5, CI.mle[2], 5, lty=lty1, lwd=lwd2, 
                      col=col1, ...)
           }
           else
             segments(CI.mle[1], 4, CI.mle[2], 4, lty=lty1, lwd=lwd2, 
                      col=col1, ...)
	   if(add.labs)
              title(xlab=paste("coefficient of", offset),
                    ylab="confidence interval",
                    main=paste(100*(1-alpha), 
                               "% Confidence Intervals "),
                    ...)
           if(add.leg == TRUE)
           {
             text(coeff[2,1], 1.2,
                  "modified likelihood root (cont.corr.)", cex=cex, ...)
             text(coeff[2,1], 2.2, "modified likelihood root",
                  cex=cex, ...)
             text(coeff[2,1], 3.2, "likelihood root", cex=cex, ...)
             if( !is.scalar )
             {
               text(coeff[2,1], 4.2, "Wald pivot (cond .MLE)", 
                    cex=cex, ...)
               text(coeff[2,1], 5.2, "Wald pivot",
                     cex=cex, ...)
             }
             else
               text(coeff[2,1], 4.2, "Wald pivot",
                     cex=cex, ...)
           }
##               ---------------
                },
                "8"={
                        par(pty="s")
                        if(is.scalar) split.screen(figs=c(1,2))
                        else split.screen(figs=c(2,2))
##               ---------------
           screen(1)
           plot(inf, type="l", xlab="", ylab="", lty=lty1, 
                lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, col=col1, ...)
           points(limINF, pch="o", cex=cex.axis, col=col1, ...)
           abline(h=0.2, lty="dotted", ...)
           abline(h=-0.2, lty="dotted", ...)
           if(add.labs)
             title(xlab=paste("coefficient of", offset), ylab="INF",
                   main="INF correction term", ...)
           screen(2)
           plot(inf.rp, type="l", xlab="", ylab="", lty=lty1, 
                   lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                cex.main=cex.main, col=col1, ...)
           if(add.labs)
             title(xlab="likelihood root", ylab="INF",
                   main="INF correction term", ...)
           if( !is.scalar )
           {
             screen(3)
             plot(np, type ="l", xlab="", ylab="", lty=lty1, 
                  lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                  cex.main=cex.main, col=col1, ...)
             points(limNP, pch="o", cex=cex.axis, col=col1, ...)
             abline(h=0.2, lty=3, ...)
             abline(h=-0.2, lty=3, ...)
             if(add.labs)
               title(xlab=paste("coefficient of", offset), ylab="NP",
                     main="NP correction term", ...)
             screen(4)
             plot(np.rp, type="l", xlab="", ylab="", lty=lty1, 
                  lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                  cex.main=cex.main, col=col1, ...)
             if(add.labs)
               title(xlab="likelihood root", ylab="NP",
                     main="NP correction term", ...)
           }
##               ---------------
                close.screen(all.screens=TRUE)
                par(pty="m")
                })
    if( missing(which) )
      pick <- menu(tmenu, 
                   title="\n Make a plot selection (or 0 to exit)\n")
    if( (pick == 0) || !missing(which) )
    {
      invisible(close.screen(all.screens=TRUE))
      break
    }
  }
}

family.cond <- function(object, ...) object$family

family.summary.cond <- function(object, ...) object$family

print.cond <- function(x, digits = max(3, getOption("digits")-3), ...)
{
  is.scalar <- x$is.scalar
  if( !is.null(cl <- x$call) )
  {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nFormula:  ")
  dput(x$formula)
  cat("Family:  ")
  print(x$family)
  cat("Offset:  ")
  if( x$offset != "1" )
    print(as.name(x$offset))
  else cat("Intercept\n")
  cat("\n")
  if( is.scalar )
    print(x$coef[1,], digits=digits)
  else print(x$coef, digits=digits)
  xd <- x$diagnostics
  cat("\nDiagnostics: \n")
  if(is.scalar)
    print(xd[1], digits=digits)
  else print(xd, digits=digits)
  cat("\n Approximation based on", x$n.approx, "points\n")
}

print.summary.cond <- function(x, all = x$all, Coef = x$cf, 
                               int = x$int, test = x$hyp,
                               digits = if(!is.null(x$digits)) x$digits
                                        else max(3, getOption("digits")-3),
                               ...)
{
  cat("\n Formula:  ")
  dput(x$formula)
  cat(" Family:  ")
  print(x$family)
  cat(" Offset:  ")
  if( x$offset != "1" )
    print(as.name(x$offset))
  else cat("Intercept\n")
  cat("\n")
  is.scalar <- x$is.scalar
  names.stat <- c("Wald pivot                             ",
                  "Wald pivot (cond. MLE)                 ",
                  "Likelihood root                        ",
                  "Modified likelihood root               ",
                  "Modified likelihood root (cont. corr.) ",
                if(all)
                  c("Lugannani-Rice approx.               ",
                    "Lugannani-Rice approx. (cont. corr.) "))
  if(Coef)
  {
    if( !is.scalar )
      print(signif(x$coefficients, digits=digits))
    else print(signif(x$coefficients[1,], digits=digits))
  }
  if(int)
  {
    dim.alpha <- length(x$alpha)
    cat("\nConfidence intervals",
        if(all) c("and one-sided confidence bounds"))
    cat("\n--------------------",
        if(all) c("-------------------------------"), sep="-")
    for(i in seq(along = x$alpha))
    {
      cat("\n level =", 100*(1-x$alpha[i]), "%\n")
      if(all)
      {
        n <- length(names.stat)
        idx <- 1:n
        ntimes <- n
        if( is.scalar )
        {
          idx <- idx[-2]
          ntimes <- ntimes - 1
        }
        print.mat <- data.frame(
                as.vector(signif(x$conf.int[i, idx],
                                 digits=digits)), rep("", ntimes),
                as.vector(signif(x$conf.int[i+2*dim.alpha, idx],
                                 digits=digits)))
        dimnames(print.mat)[[2]] <- c("lower", "two-sided", "upper")
        dimnames(print.mat)[[1]] <- names.stat[idx]
        print(print.mat)
        cat("\n")
        print.mat <- data.frame(
                as.vector(signif(x$conf.int[dim.alpha+i, idx],  
                                digits=digits)), rep("", ntimes),
                as.vector(signif(x$conf.int[i+3*dim.alpha, idx],
                                 digits=digits)))
        dimnames(print.mat)[[2]] <- c("lower", " one-sided", "upper")
        dimnames(print.mat)[[1]] <- names.stat[idx]
        print(print.mat)
      }
      else
      {
        n <- length(names.stat)
        idx <- 1:n
        ntimes <- n
        if(is.scalar)
        {
          idx <- idx[-2]
          ntimes <- ntimes - 1
        }
        print.mat <- data.frame(
                as.vector(signif(x$conf.int[i, idx],
                                 digits=digits)), rep("", ntimes),
                as.vector(signif(x$conf.int[i+2*dim.alpha, idx],
                                 digits=digits)))
        dimnames(print.mat)[[2]] <- c("lower", "two-sided", "upper")
        dimnames(print.mat)[[1]] <- names.stat[idx]
        print(print.mat)
      }
    }
  }
  qT <- x$signif.tests$qTerm
  if( !is.null(test) )
  {
    dim.test <- length(x$hypotheses)
    cat("\nTest statistics")
    cat("\n---------------")
    for(i in seq(along = x$hypotheses))
    {
      mat <- matrix(signif(x$signif.tests$stats[1:(ifelse(all, 7, 5)),
                                                c(i, i+dim.test)],
                           digits = digits), ncol = 2)
      dimnames(mat) <- list(names.stat, c(" statistic ",
                                          " tail prob. "))
      txt <- if(x$offset == "1")
               "\n hypothesis : Intercept ="
             else paste("\n hypothesis : coef(", x$offset, ") =")
      txt <- paste(txt, x$hypotheses[i], "\n")
      cat(txt)
      if( !is.scalar )
        print(mat)
      else print(mat[-2,])
      xq <- matrix(qT[i,], nrow=1)
      dimnames(xq) <- list("\"q\" correction term:", "")
      print(xq, digits=digits)
    }
  }
  xd <- x$diagnostics
  cat("\nDiagnostics:")
  cat("\n----------- \n")
  if( !is.scalar )
    print(xd, digits=digits)
  else print(xd[1], digits=digits)
  cat("\n Approximation based on", x$n.approx, "points\n")
}
