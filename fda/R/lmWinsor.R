lmWinsor <- function(formula, data, lower=NULL, upper=NULL, trim=0,
        quantileType=7, subset, weights=NULL, na.action,
        model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
        singular.ok = TRUE, contrasts = NULL, offset=NULL,
        method=c('QP', 'clip'), eps=sqrt(.Machine$double.eps),
        trace=ceiling(sqrt(nrow(data))), ...)
{
##
## 1.  Identify inputs and outputs
##
  cl <- match.call()
  if(missing(na.action))
    na.action <- get(options('na.action')$na.action)
  mdly <- mdlx <- formula
  mdly[[3]] <- NULL
  mdlx[[2]] <- NULL
  xNames <- all.vars(mdlx)
  yName <- all.vars(mdly)
  if(length(yName) != 1){
    mdl. <- paste(as.character(formula)[c(2, 1, 3)], collapse=" ")
    msg <- paste("'formula' must include a single response;  found",
                 length(yName), 'in', mdl.)
    stop(msg)
  }
  if(as.character(mdly[[2]]) != yName)
    stop("lmWinsor can not accept a formula with a transformed 'y';",
         "  left hand side = ", mdly[[2]], ";  y = ", yName)
##
## 2.  output = a modified 'lm' object or a list of such objects?
##
  N <- nrow(data)
  if(missing(subset))subset <- 1:N
  if(is.null(weights))weights <- rep(1, N)
#
  if(is.list(lower) || is.list(upper) || (length(trim)>1)){
    fit <- lmWinsor2(formula=formula, data=data, lower=lower,
        upper=upper, trim=trim, quantileType=quantileType, subset=subset,
        weights=weights, na.action=na.action, model = model, x = x, y = y,
        qr = qr, singular.ok = singular.ok, contrasts = contrasts,
        offset=offset, method=method, eps=eps, trace=trace, ... )
    fit$call <- cl
    return(fit)
  }
#
  fit <- lmWinsor1(formula=formula, data=data, lower=lower,
        upper=upper, trim=trim, quantileType=quantileType, subset=subset,
        weights=weights, na.action=na.action, model = model, x = x, y = y,
        qr = qr, singular.ok = singular.ok, contrasts = contrasts,
        offset=offset, method=method, eps=eps, trace=trace, ...)
  fit$call <- cl
  fit
}

lmWinsor1 <- function(formula, data, lower=NULL, upper=NULL, trim=0,
        quantileType=7, subset, weights=NULL, na.action,
        model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
        singular.ok = TRUE, contrasts = NULL, offset=NULL,
        method=c('QP', 'clip'), eps=sqrt(.Machine$double.eps),
        trace=ceiling(sqrt(nrow(data))), ...){
##
## 1.  yName, xNames
##
  start.time <- proc.time()
  cl <- match.call()
  mdly <- mdlx <- formula
  mdly[[3]] <- NULL
  mdlx[[2]] <- NULL
  xNames <- all.vars(mdlx)
  yName <- all.vars(mdly)
##
## 2.  Check 'lower' and 'upper'
##
#  2.1.  Do lower and upper have names?
  lowNames <- names(lower)
  if(length(lowNames) != length(lower))
    stop("lower must have names")
  hiNames <- names(upper)
  if(length(hiNames) != length(upper))
    stop("upper must have names")
#
  goodl <- (lowNames %in% names(data))
  if(!all(goodl))
    warning("names(lower) not in names(data);  first one = ",
            lowNames[goodl][1])
  goodh <- (hiNames %in% names(data))
  if(!all(goodh))
    warning("names(upper) not in names(data);  first one = ",
            hiNames[goodh][1])
#  2.2.  Identify numeric columns of 'data'
  clipData <- na.action(data[unique(c(yName, xNames))])
#
  numVars <- sapply(clipData, is.numeric)
  numV <- names(numVars[numVars])
#  2.3.  Are numeric variables in lower and upper?
  numLower <- (numV %in% names(lower))
  numUpper <- (numV %in% names(upper))
  nnV <- length(numV)
#  2.4.  Some numeric variables are not in lower;  add
  {
    if(!all(numLower)){
      Lower <- rep(NA, nnV)
      names(Lower) <- numV
      gotLo <- (numV %in% lowNames)
      if(any(gotLo)) {
        loGot <- lower[numV[gotLo]]
        Lower[gotLo] <- loGot
      }
      for(v in numV[!numLower])
        Lower[v] <- quantile(clipData[[v]], trim, na.rm=TRUE, names=FALSE,
                             type=quantileType)
    }
    else
      Lower <- lower
  }
  {
    if(!all(numUpper)){
      Upper <- rep(NA, nnV)
      names(Upper) <- numV
      gotHi <- (numV %in% hiNames)
      if(any(gotHi)){
        hiGot <- upper[numV[gotHi]]
        Upper[gotHi] <- hiGot
      }
      for(v in numV[!numUpper])
        Upper[v] <- quantile(clipData[[v]], 1-trim, na.rm=TRUE, names=FALSE,
                             type=quantileType)
    }
    else
      Upper <- upper
  }
##
## 3.  clipData = data with xNames clipped to (Lower, Upper)
##
  for(x. in intersect(xNames, numV)){
    x.L <- Lower[x.]
    xl <- pmax(data[[x.]], x.L)
#    xl <- pmax(data[[x.]], x.L*(1+3*.Machine$double.eps))
    x.U <- Upper[x.]
    clipData[[x.]] <- pmin(xl, x.U)
#    clipData[[x.]] <- pmin(xl, x.U *(1-3*.Machine$double.neg.esp))
  }
##
## 4.  fit <- lm(...)
##
  N <- nrow(clipData)
# if(missing(subset))subset <- 1:N
#
  cl0 <- as.list(cl)
  cl0[[1]] <- NULL
  cl0$lower <- NULL
  cl0$upper <- NULL
  cl0$trim <- NULL
  cl0$quantileType <- NULL
  cl0$eps <- NULL
  cl0$trac <- NULL
  cl0$method <- NULL
  cl0$trace <- NULL
  if(is.null(offset))cl0$offset <- NULL
  cl0$data <- as.name('clipData')
  cl0$formula <- eval(cl0$formula)
  cl0$subset <- eval(cl0$subset)
  cl0$weights <- eval(cl0$weights)
#
  fit <- do.call('lm', cl0)
##
## 5.  Convert to class 'lmWinsor'
##
  fit$call <- cl
  fit$lower <- Lower
  fit$upper <- Upper
  pred <- predict.lm(fit)
  class(fit) <- c("lmWinsor", class(fit))
##
## 6.  all fit$fitted.values %in% (Lower, Upper)[yName]?
##
## Need yL < yL. < yLin <= yUin < yU. < yU
## with (yLin-yL)~=2*eps & (yU-yUin)~=2*eps
##
  y. <- data[[yName]]
  mod <- mean(abs(y.[abs(y.)<Inf]))
  Eps <- {
    if(mod==0) eps
    else eps*mod
  }
#  yL <- (Lower[yName]*(1-3*.Machine$double.neg.eps))
  yLin <- Lower[yName]+Eps
#  yU <- (Upper[yName]*(1+3*.Machine$double.eps))
  yUin <- Upper[yName]-Eps
# yL == yU?
  if(yLin>yUin){
    yLin <- yUin <- (yLin+yUin)/2
  }
  yL. <- yLin-Eps
  yL <- yL.-Eps
  yLow <- (pred<yL)
#
  yU. <- yUin+Eps
  yU <- yU.+Eps
  yHi <- (yU<pred)
#
  mtd <- match.arg(method)
  fit$out <- cbind(low=yLow, high=yHi)
  if((mtd=='clip') || !any(yLow | yHi)){
    fit$fitted.values <- pmax(yL., pmin(yU., pred))
    if(!model) fit$model <- NULL
    fit$message <- {
      if(mtd=='clip')
        paste(sum(yLow | yHi), "of", N, "fitted.values clipped")
      else 'Initial fit in bounds'
    }
    fit$elapsed.time <- max(proc.time()-start.time, na.rm=TRUE)
    return(fit)
  }
##
## 7.  Else use quadratic progamming to minimize the
##     Winsorized sum of squares of residuals
##
  if(!require(quadprog)){
    msg <- paste(sum(yLow | yHi), "of", nrow(data),
                 "predictions outside bounds, but package quadprog",
                 "not installed;\ncan not iterate;  using method='clip'")
    warning(msg)
    fit$fitted.values <- pmax(yL., pmin(yU., pred))
    if(!model) fit$model <- NULL
    fit$message <- msg
    return(fit)
  }
#
  out <- matrix(FALSE, N, 2, dimnames=list(NULL,
                       c('out', 'high')))
  coef0 <- coef(fit)
  k <- sum(!is.na(coef0))
  extraStats <- c("newConstraint", "SSEraw", "SSEclipped",
      "nLoOut", "nLo.", "nIn", "nHi.", "nHiOut")
  ks <- length(extraStats)
  coefiter <- matrix(NA, N+1, k+ks, dimnames=list(NULL,
        c(names(coef0[!is.na(coef0)]), extraStats)) )
  coefiter[1, 1:k] <- coef0[!is.na(coef0)]
  coefiter[1, "newConstraint"] <- 0
  coefiter[1, "SSEraw"] <- sum((y.-pred)^2)
  predW <- pmax(yL., pmin(yU., pred))
  coefiter[1, "SSEclipped"] <- sum((y.-predW)^2)
  coefiter[1, "nLoOut"] <- sum(yLow)
  coefiter[1, "nLo."] <- sum(!yLow & (pred<yLin))
  coefiter[1, "nHiOut"] <- sum(yHi)
  coefiter[1, "nHi."] <- sum(!yHi & (yUin<=pred))
  coefiter[1, "nIn"] <- (N-sum(coefiter[1,
           c("nLoOut", "nLo.", "nHi.", "nHiOut") ]))
#
#  templimits <- matrix(NA, N, 3, dimnames=list(NULL,
#                       c("newConstraint", "lower", "upper") ) )
#
  X <- model.matrix(fit)[, !is.na(coef0), drop=FALSE]
#  yMin <- min(y.)
#  yMax <- max(y.)
#
  w.5 <- sqrt(weights)
  for(i in 1:N){
#   The standard exit from this loop is via 'break'
#   afer all constraints are satisfied.
#   7.1.  Find prediction farthest out
    dLow <- (yL-pred)
    dHi <- (pred-yU)
    dLow. <- dLow[!out[, 'out']]
    dHi. <- dHi[!out[, 'out']]
    i1 <- i+1
    {
      if(max(dLow.) > max(dHi.)){
        lows <- (which(!out[, 'out'])[max(dLow.)==dLow.])
        lowest <- lows[y.[lows]==min(y.[lows])][1]
        out[lowest, 'out'] <- TRUE
#        templimits[i, "newConstraint"] <- lowest
        coefiter[i1, "newConstraint"] <- lowest
#        yMin <- min(yLin, pred[pred>pred[lowest]])
      }
      else{
        highs <- (which(!out[, 'out'])[max(dHi.)==dHi.])
        highest <- highs[y.[highs]==max(y.[highs])][1]
        out[highest, ] <- TRUE
#       templimits[i, "newConstraint"] <- highest
        coefiter[i1, "newConstraint"] <- highest
#        yMax <- max(yUin, pred[pred<pred[highest]])
      }
    }
#   7.2.  Use QP ...
#    Dmat. <- crossprod(X[!out[, 1], , drop=FALSE])
    in. <- !out[, 'out']
    if(sum(in.)<1){
      fit$message <- paste('Quadratic programming iteration ended',
                           'with all predictions outside limits.')
      fit$qr <- qrX.old
      break
    }
    qrX <- qr(w.5[in.]*X[in., , drop=FALSE])
    qR <- qr.R(qrX)
    signR <- sign(diag(qR))
    qR. <- (signR * qR)
#    qrXdiag <- diag(qR.)
    qrXrng <- range(diag(qR.))
    sing <- ((dim(qR)[1] < k) || (qrXrng[1] < eps * qrXrng[2]))
    if(sing){
      fit$message <- 'Iteration terminated by a singular quadratic program'
      fit$qr <- {
        if(i<2) qrX else qrX.old
      }
      break
    }
    Dmat. <- solve(qR.)
#
    dvec. <- as.numeric(crossprod(X[in., , drop=FALSE],
                                  weights[in.]*y.[in.]))
#
    outL <- (out[, 'out'] & !out[, 'high'])
    outH <- (out[, 'out'] & out[, 'high'])
    AmatL <- X[outL, , drop=FALSE]
    AmatH <- X[outH, , drop=FALSE]
    Amat. <- rbind(-AmatL, AmatH)
#
#    maxPredLo <- max(pred[outL], -Inf)
#    minPredHi <- min(pred[outH], Inf)
#    yMin <- min(y.[y.>maxPredLo], yL)+Eps
#
#    yMax <- max(y.[y.<minPredHi], yU)-Eps
#
#    bvec. <- c(rep(-yMin, nrow(AmatL)), rep(yMax, nrow(AmatH)))
    bvec. <- c(rep(-yL., nrow(AmatL)), rep(yU., nrow(AmatH)))
    Ab <- constraintCheck(t(Amat.), bvec.)
#
#    templimits[i, c("lower", "upper")] <- c(yMin, yMax)
#    templimits[i, c("lower", "upper")] <- c(yL., yU.)
#
#    QPi <- solve.QP(Dmat=Dmat., dvec=dvec., Amat=Amat., bvec=bvec.)
    env <- new.env()
    assign("D", Dmat., envir=env)
    assign("d", dvec., envir=env)
#    assign("A", t(Amat.), envir=env)
#    assign("b", bvec., envir=env)
    assign("A", Ab$A, envir=env)
    assign("b", Ab$b, envir=env)
    QPlist <- list(Dmat=quote(D), dvec=quote(d),
                   Amat=quote(A), bvec=quote(b),
                   factorized=TRUE)
    QPi <- do.call('solve.QP', QPlist, envir=env)
#    QPi <- solve.QP(Dmat., as.numeric(dvec), Amat, bvec)
#   7.3.  Are unconstrained predictions inside?
#    pred.old <- pred
#    i1 <- i+1
    coefiter[i1, 1:k] <- QPi$solution
    pred <- X %*% QPi$solution
#    pred <- pmax(yL, pmin(yU, pred0))
    yLow <- (pred<yL)
    yHi <- (yU<pred)
    coefiter[i1, "SSEraw"] <- sum((y.-pred)^2)
    predW <- pmax(yL., pmin(yU., pred))
    coefiter[i1, "SSEclipped"] <- sum((y.-predW)^2)
    coefiter[i1, "nLoOut"] <- sum(yLow)
    coefiter[i1, "nLo."] <- sum(!yLow & (pred<yLin))
    coefiter[i1, "nHiOut"] <- sum(yHi)
    coefiter[i1, "nHi."] <- sum(!yHi & (yUin<=pred))
    coefiter[i1, "nIn"] <- (N-sum(coefiter[i1,
           c("nLoOut", "nLo.", "nHi.", "nHiOut") ]))
#
    out2 <- (yLow | yHi)
    if(!any(out2[!out[, 'out']])) {
      fit$message <- 'QP iterations successful'
      fit$qr <- qrX
      break
    }
    qrX.old <- qrX
    if((trace > 0) && ((i%%trace) == 0))cat(i, "")
#   7.4.  Find the next extreme prediction ... -> 7.1
  }
  if((trace>0) && (i >= trace))cat('\n')
##
## 8.  Modify fit as appropriate
##
  coef0[!is.na(coef0)] <- QPi$solution
  fit$coefficients <- coef0
#
  predy <- pmax(yL, pmin(yU, pred))
  fit$residuals <- (y.-predy)
  fit$fitted.values <- predy
  fit$coefIter <- coefiter[!is.na(coefiter[, 1]), , drop=FALSE]
#  fit$tempLimits <- templimits[!is.na(templimits[, 1]), , drop=FALSE]
#
  fit$elapsed.time <- max(proc.time()-start.time, na.rm=TRUE)
  fit
}

constraintCheck <- function(Amat, bvec,
           eps=sqrt(.Machine$double.eps) ){
##
## 1.  Make cols of Amat^2 sum to 1
##
  a. <- sqrt(apply(Amat^2, 2, sum))
  n <- nrow(Amat)
  a <- (Amat / rep(a., each=n))
  b. <- (bvec / a.)
##
## 2.  Check for identical rows in 'a'
##
  d.a <- dist(t(a))
  if(any(d.a<eps)){
    d.a2 <- as.matrix(d.a)
    diag(d.a2) <- 1
    rr <- row(d.a2)[d.a2<eps][1]
    cr <- col(d.a2)[d.a2<eps][1]
    {
      if(b.[rr]<b.[cr])
        return(constraintCheck(Amat[, -rr, drop=FALSE], bvec[-rr]))
      else
        return(constraintCheck(Amat[, -cr, drop=FALSE], bvec[-cr]))
    }
  }
##
## 3.  No redundancies found
##
  list(A=Amat, b=bvec)
}

lmWinsor2 <- function(formula, data, lower=NULL, upper=NULL, trim=0,
        quantileType=7, subset, weights=NULL, na.action,
        model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
        singular.ok = TRUE, contrasts = NULL, offset=NULL,
        method=c('QP', 'clip'), eps=sqrt(.Machine$double.eps),
        trace=ceiling(sqrt(nrow(data))), ...){
##
## 1.  Rationalize lower, upper, trim
##
  nt <- length(trim)
#
  {
    if(is.null(lower))
      nl <- 0
    else {
      if(!is.list(lower))
        lower <- list(lower)
      nl <- length(lower)
    }
  }
#
  {
    if(is.null(upper))
      nu <- 0
    else {
      if(!is.list(upper))
        upper <- list(upper)
      nu <- length(upper)
    }
  }
#
  n. <- max(nt, nl, nu)
  trim <- rep(trim, length=n.)
  if((0 < nl) & (nl < n.))
    for(i in 1:(n.-nl))
      lower[[nl+i]] <- lower[[i]]
#
  if((0 < nu) & (nu < n.))
    for(i in 1:(n.-nu))
      upper[[nu+i]] <- upper[[i]]
# preserve names if any
  limNames <- NULL
  if(nt == n.)
    limNames <- names(trim)
  if(is.null(limNames) && (nl == n.))
    limNames <- names(lower)
  if(is.null(limNames) && (nu == n.))
    limNames <- names(upper)
##
## 2.  Do it
##
  lmW <- vector('list', n.)
  names(lmW) <- limNames
#
  for(i in 1:n.)
    lmW[[i]] <- lmWinsor1(formula, data, lower=lower[[i]],
        upper=upper[[i]], trim=trim[i], quantileType=quantileType,
        subset=subset, weights=weights, na.action=na.action,
        model = model, x = x, y = y, qr = qr,
        singular.ok = singular.ok, contrasts = contrasts,
        offset=offset, method=method, eps=eps, trace=trace, ...)
##
## 3.  Done
##
  class(lmW) <- "lmWinsor"
  lmW
}
