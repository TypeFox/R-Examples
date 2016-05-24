globalVariables(".SD")

logit <- function(x) log(x) - log(1 - x)
odds <- function(x) x / (1 - x)
  
is.binary <- function(v){

  if(is.numeric(v) & all(v == 0 | v == 1, na.rm=TRUE))
    TRUE
  else
    FALSE
    
}


stdGlm <- function(fit, data, X, x, clusters, case.control=FALSE){

  #---PREPARATION---

  formula <- fit$formula
  weights <- fit$prior.weights
  npar <- length(fit$coef)
  out <- list(fit=fit, X=X)

  #delete rows with missing on variables in the model
  rownames(data) <- 1:nrow(data)
  m <- model.frame(formula=formula, data=data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  n <- nrow(data)

  #Can write code more generally with
  #if(missing(clusters)) clusters <- 1:nrow(data)
  #but 2 problems when constructing meat in sandwich formula:
  #1) must always aggregate, which is a bit slow, even though much faster
  #when using data.table than the aggregate function, 2) still need
  #to have special case for non-clustered case-control

  if(!missing(clusters)){
    data <- data[order(data[, clusters]), ]
    clusters <- data[, clusters]
    n.cluster <- length(unique(clusters))
  }

  #assign values to x and reference if not supplied
  #make sure x is a factor if data[, X] is a factor
  if(missing(x)){
    if(is.factor(data[, X]))
      x <- as.factor(levels(data[, X]))
    if(is.numeric(data[, X]))
      if(is.binary(data[, X]))
        x <- c(0, 1)
      else
        x <- round(mean(data[, X], na.rm=TRUE), 2)
  }
  else{
    if(is.factor(x)){
      temp <- x
      levels(x) <- levels(data[, X])
      x[1:length(x)] <- temp
    }
    else{
      if(is.factor(data[, X])){
        x <- factor(x)
        temp <- x
        levels(x) <- levels(data[, X])
        x[1:length(x)] <- temp
      }
      else 
        x <- sort(x)
    }
  }
  
  nX <- length(x)
  out <- c(out, list(x=x))

  #---ESTIMATES OF MEANS AT VALUES SPECIFIED BY x ---
  
  pred <- matrix(nrow=n, ncol=nX)
  for(i in 1:nX){
    data.temp <- data
    data.temp[, X] <- x[i]
    pred[, i] <- predict(object=fit, newdata=data.temp, type="response")
  }
  means.est <- colSums(weights * pred, na.rm=TRUE) / sum(weights)
  
  #---VARIANCE OF MEANS AT VALUES SPECIFIED BY x---
 
  ores <- weights * residuals(object=fit, type="response") *
    model.matrix(object=formula, data=data)
  mres <- weights * (pred - matrix(rep(means.est, each=n), nrow=n, ncol=nX))
  res <- cbind(mres, ores)
  
  if(case.control & missing(clusters)){
    outcome <- as.character(formula)[2]
    controls <- which(data[, outcome]==0)
    cases <- which(data[, outcome]==1)
    n0 <- length(controls)
    n1 <- length(cases)
    J <- n0/n * var(res[controls, ], na.rm=TRUE) + n1/n * var(res[cases, ], na.rm=TRUE)
  }
  else{
    if(!missing(clusters)){
      res <- data.table(res)
      res <- res[, j=lapply(.SD,sum), by=clusters]
      res <- as.matrix(res)
      res <- res[, -1]
    }
    J <- var(res, na.rm=TRUE)
  }
  
  g <- family(fit)$mu.eta
  SI.beta <- matrix(nrow=nX, ncol=npar)
  for(i in 1:nX){
    data.temp <- data
    data.temp[, X] <- x[i]
    dmu.deta <- g(predict(object=fit, newdata=data.temp))
    deta.dbeta <- model.matrix(object=delete.response(terms(fit)), data=data.temp)
    dmu.dbeta <- dmu.deta * deta.dbeta
    SI.beta[i, ] <- colMeans(weights * dmu.dbeta)
  }
  SI <- cbind(-diag(nX) * mean(weights), SI.beta)
  oI <- cbind(matrix(0, nrow=npar, ncol=nX), -solve(vcov(object=fit)) / n)
  I <- rbind(SI, oI)
  if(missing(clusters))
    V <- (solve(I) %*% J %*% t(solve(I)) / n)[1:nX, 1:nX]
  else
    V <- (solve(I) %*% J %*% t(solve(I)) * n.cluster / n^2)[1:nX, 1:nX]
  means.vcov <- V

  out <- c(out, list(means.est=means.est, means.vcov=means.vcov))

  #---OUTPUT---

  class(out) <- "stdGlm"
  return(out)

}

summary.stdGlm <- function(object, CI.type="plain", CI.level=95,
  transform=NULL, contrast=NULL, reference=NULL, ...){

  qqq <- qnorm((1+CI.level/100)/2)

  est <- object$means.est
  V <- as.matrix(object$means.vcov)
  nX <- length(object$x)

  if(!is.null(transform)){
    if(transform == "log"){
      dtransform.dm <- diag(1 / est, nrow=nX, ncol=nX)
      est <- log(est)
    }
    if(transform == "logit"){
      dtransform.dm <- diag(1 / (est * (1 - est)), nrow=nX, ncol=nX)
      est <- logit(est)
    }
    if(transform == "odds"){
      dtransform.dm <- diag(1 / (1 - est)^2, nrow=nX, ncol=nX)
      est <- odds(est)
    }
    V <- t(dtransform.dm) %*% V %*% dtransform.dm
  }

  if(!is.null(contrast)){
    referencepos <- match(reference, object$x)
    if(contrast == "difference"){
      dcontrast.dtransform <- diag(nX)
      dcontrast.dtransform[referencepos, ] <- -1
      dcontrast.dtransform[referencepos, referencepos] <- 0
      est <- est-est[referencepos]
    }
    if(contrast == "ratio"){
      dcontrast.dtransform <- diag(1 / est[referencepos], nrow=nX, ncol=nX)
      dcontrast.dtransform[referencepos, ] <- -est / est[referencepos]^2
      dcontrast.dtransform[referencepos, referencepos] <- 1
      est <- est / est[referencepos]
    }
    V <- t(dcontrast.dtransform) %*% V %*% dcontrast.dtransform
    V[referencepos, ] <- 0
    V[, referencepos] <- 0
  }

  se <-  sqrt(diag(V))

  if(CI.type == "plain"){
    theta <- qqq * se
    lower <- est - theta
    upper <- est + theta
  }
  if(CI.type == "log"){
      theta <- qqq * se / est
      lower <- est * exp(-theta)
      upper <- est * exp(theta)
  }
  if(is.factor(reference))
    reference <- as.character(reference)
  est.table <- as.matrix(cbind(est, se, lower, upper), nrow=length(est), ncol=4)
  dimnames(est.table) <- list(object$x,
    c("Estimate", "Std. Error", paste("lower",CI.level), paste("upper",CI.level)))
  out <- c(object, list(est.table=est.table,transform=transform,contrast=contrast,reference=reference))

  class(out) <- "summary.stdGlm"
  return(out)

}

print.summary.stdGlm <- function(x, ...){

  cat("\nFormula: ")
  print(x$fit$formula)
  cat("Family:",  x$fit$family$family,  "\n")
  cat("Link function:",  x$fit$family$link,  "\n")
  cat("Exposure: ", x$X,  "\n")
  if(!is.null(x$transform))
    cat("Transform: ", x$transform,  "\n")
  if(!is.null(x$contrast)){
    cat("Reference level: ", x$X, "=", x$reference,  "\n")
    cat("Contrast: ", x$contrast,  "\n")
  }
  cat("\n")
  print(x$est.table, digits=3)

}

plot.stdGlm <- function(x, CI.type="plain", CI.level=95,
  transform=NULL, contrast=NULL, reference=NULL, ...){

  object <- x
  x <- object$x

  dots <- list(...)

  xlab <- object$X
  
  if(is.factor(reference))
    reference <- as.character(reference)

  if(is.null(contrast)){
    if(is.null(transform))
      ylab <- expression(mu)
    else{
      if(transform == "log")
        ylab <- expression(paste("log(", mu, ")"))
      if(transform == "logit")
        ylab <- expression(paste("logit(", mu, ")"))
      if(transform == "odds")
        ylab <- expression(paste(mu, " / (1-", mu, ")"))
    }
  }
  else{
    if(contrast == "difference"){
      if(is.null(transform))
        ylab <- c(bquote(paste(mu, " - ", mu[.(reference)])), expression())
      else{
        if(transform == "log")
          ylab <- c(bquote(paste(log, "(", mu, ") - ", log, "(",
            mu[.(reference)], ")", sep="")), expression())
        if(transform == "logit")
          ylab <- c(bquote(paste(logit, "(", mu, ") - ", logit,
           "(", mu[.(reference)], ")", sep="")), expression())
        if(transform == "odds")
          ylab <- c(bquote(paste(mu, " / (", 1 - mu, ") - ",
            mu[.(reference)], " / (", 1 - mu[.(reference)], ")", sep="")),
            expression())
      }
    }
    if(contrast == "ratio"){
      if(is.null(transform))
        ylab <- c(bquote(paste(mu, " / ", mu[.(reference)])), expression())
      else{
        if(transform == "log")
          ylab <- c(bquote(paste(log, "(", mu, ") / ", log, "(",
            mu[.(reference)], ")", sep="")), expression())
        if(transform == "logit")
          ylab <- c(bquote(paste(logit, "(", mu, ") / ", logit,
            "(", mu[.(reference)], ")", sep="")), expression())
        if(transform == "odds")
          ylab <- c(bquote(paste(mu, " / (", 1 - mu, ") / ",
            mu[.(reference)], " / (", 1 - mu[.(reference)],
            ")", sep="")), expression())
      }
    }
  }

  sum.obj <- summary(object=object, CI.type=CI.type, CI.level=CI.level,
    transform=transform, contrast=contrast, reference=reference)
  est <- sum.obj$est.table[, 1]
  lower <- sum.obj$est.table[, 3]
  upper <- sum.obj$est.table[, 4]
  
  ylim <- c(min(c(lower,upper)), max(c(lower,upper)))
  
  if(is.numeric(x) & length(x)>1){
    args <- list(x=x, y=x, xlab=xlab, ylab=ylab, ylim=ylim, type="n")
    args[names(dots)] <- dots
    do.call("plot", args=args)
    lines(x, est)
    lines(x, upper, lty=3)
    lines(x, lower, lty=3)
  }
  if(is.factor(x) | is.binary(x) | (is.numeric(x) & length(x)==1)){
    args <- list(x=1:length(x), y=1:length(x), xlab=xlab, ylab=ylab,
      xlim=c(0, length(x)+1), ylim=ylim, type="n", xaxt="n")
    args[names(dots)] <- dots
    do.call("plot", args=args)
    points(1:length(x), est)
    points(1:length(x), upper, pch=0)
    points(1:length(x), lower, pch=0)
    for(i in 1:length(x))
      lines(x=c(i, i), y=c(lower[i], upper[i]), lty="dashed")
    mtext(text=x, side=1, at=1:length(x))
  }
}

stdCoxph <- function(fit, data, X, x, t, clusters){

  #---PREPARATION---

  formula <- fit$formula
  npar <- length(fit$coef)
  fit.detail <- coxph.detail(object=fit)
  out <- list(fit=fit, X=X)

  #delete rows with missing on variables in the model 
  rownames(data) <- 1:nrow(data)
  m <- model.frame(formula=formula, data=data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  n <- nrow(data)
  
  #extract end variable and event variable
  Y <- model.extract(frame=m, "response")
  if(ncol(Y) == 2){
    end <- Y[, 1]
    event <- Y[, 2]
  }
  if(ncol(Y) == 3){
    end <- Y[, 2]
    event <- Y[, 3]
  }
  
  if(is.null(fit$weights))
    weights <- rep(1, nrow(data))
  else
    weights <- fit$weights

  #Can write code more generally with
  #if(missing(clusters)) clusters <- 1:nrow(data)
  #but a problem when constructing meat in sandwich formula:
  #must always aggregate, which is a bit slow, even though much faster
  #when using data.table than the aggregate function

  if(!missing(clusters)){
    data <- data[order(data[, clusters]), ]
    clusters <- data[, clusters]
    n.cluster <- length(unique(clusters))
  }
  
  #assign values to x and reference if not supplied
  #make sure x is a factor if data[, X] is a factor
  if(missing(x)){
    if(is.factor(data[, X]))
      x <- as.factor(levels(data[, X]))
    if(is.numeric(data[, X]))
      if(is.binary(data[, X]))
        x <- c(0, 1)
      else
        x <- round(mean(data[, X], na.rm=TRUE), 2)
  }
  else{
    if(is.factor(x)){
      temp <- x
      levels(x) <- levels(data[, X])
      x[1:length(x)] <- temp
    }
    else{
      if(is.factor(data[, X])){
        x <- factor(x)
        temp <- x
        levels(x) <- levels(data[, X])
        x[1:length(x)] <- temp
      }
      else 
        x <- sort(x)
    }
  }
  
  nX <- length(x)    
  out <- c(out, list(x=x))

  #sort on "end-variable"
  ord <- order(end)
  data <- data[ord, ]
  weights <- weights[ord]
  if(!missing(clusters))
    clusters <- clusters[ord]
  end <- end[ord]
  event <- event[ord]

  #assign value to t if missing
  if(missing(t))
    t <- fit.detail$time
  out <- c(out, list(t=t))
  nt <- length(t)

  if(sum(fit.detail$time<=min(t)) == 0)
    stop("No events before first value in t", call.=FALSE)

  #collect important stuff from copxh.detail
  surv.est <- matrix(nrow=nt, ncol=nX)
  surv.vcov <- vector(mode="list", length=nt)
  dH0 <- fit.detail$hazard
  E <- matrix(0, nrow=n, ncol=npar)
  means <- fit.detail$means
  means <- means[rep(1:nrow(means), fit.detail$nevent), ] #handle ties
  E[event == 1, ] <- means
  H0 <- cumsum(dH0)
  H0step <- stepfun(fit.detail$time, c(0, H0))
  H0i=rep(0, n)
  dH0.untied <- rep(dH0, fit.detail$nevent) / rep(fit.detail$nevent, fit.detail$nevent)
  H0i[event == 1] <- dH0.untied * n #handle ties
  betares <- weights * residuals(object=fit, type="score")
  betares <- betares[ord, ]

  #---LOOP OVER nt

  for(j in 1:nt){

    H0it <- H0i * (end <= t[j])

    #---ESTIMATES OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x, ---

    si <- matrix(nrow=n, ncol=nX)
    PredX <- matrix(nrow=n, ncol=nX)
    tempmat <- matrix(nrow=nX, ncol=npar)
    for(i in 1:nX){
      data.temp <- data
      data.temp[, X] <- x[i]
      predX <- predict(object=fit, newdata=data.temp, type="risk")
      si[, i] <- exp(-H0step(t[j]) * predX)
      PredX[, i] <- predX
      tempmat[i, ] <- colMeans(model.matrix(object=delete.response(terms(fit)), data=data.temp)[, -1] *
        predX * si[, i] * weights)
    }
    surv.est[j, ] <- colSums(weights * si, na.rm=TRUE) / sum(weights)

    #---VARIANCE OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x, ---
    
    sres <- weights * (si - matrix(rep(surv.est[j, ],each=n), nrow=n, ncol=nX))
    H0rest <-  H0it - H0step(t[j])
    res <- cbind(sres, betares, H0rest)
    if(!missing(clusters)){
      res <- data.table(res)
      res <- res[, j=lapply(.SD,sum), by=clusters]
      res <- as.matrix(res)
      res <- res[, -1]
    }
    J <- var(res, na.rm=TRUE)
    SI <- cbind(-diag(nX) * mean(weights), -tempmat * H0step(t[j]),
      -colMeans(PredX * si * weights))
    betaI <- cbind(matrix(0, nrow=npar, ncol=nX), -solve(vcov(object=fit)) / n,
      rep(0, npar))
    H0I <- c(rep(0, nX), -colMeans(E * H0it, na.rm=TRUE), -1)
    I <- rbind(SI, betaI, H0I) 
    if(missing(clusters))
      V <- (solve(I) %*% J %*% t(solve(I)) / n)[1:nX, 1:nX]
    else
      V <- (solve(I) %*% J %*% t(solve(I)) * n.cluster / n^2)[1:nX, 1:nX]
    surv.vcov[[j]] <- V

  }
  
  out <- c(out, list(surv.est=surv.est, surv.vcov=surv.vcov))

  #---OUTPUT---

  class(out) <- "stdCoxph"
  return(out)

}

summary.stdCoxph <- function(object, t, CI.type="plain", CI.level=95,
  transform=NULL, contrast=NULL, reference=NULL, ...){

  qqq <- qnorm((1+CI.level/100)/2)

  est.all <- object$surv.est
  V.all <- object$surv.vcov
  nX <- length(object$x)
  if(missing(t))
    t <- object$t
  nt <- length(t)

  est.table <- vector(mode="list", length=nt)
  for(j in 1:nt){

    if(min(abs(t[j]-object$t)) > sqrt(.Machine$double.eps))
      stop("The standardized survival function is not estimated at t", call.=FALSE)
    else
      k <- which.min(abs(t[j] - object$t))

    est <- est.all[k, ]
    V <- as.matrix(V.all[[k]])

    if(!is.null(transform)){
      if(transform == "log"){
        dtransform.dm <- diag(1 / est, nrow=nX, ncol=nX)
        est <- log(est)
      }
      if(transform == "logit"){
        dtransform.dm <- diag(1 / (est * (1 - est)), nrow=nX, ncol=nX)
        est <- logit(est)
      }
      if(transform == "odds"){
        dtransform.dm <- diag(1 / (1 - est)^2, nrow=nX, ncol=nX)
        est <- odds(est)
      }
      V <- t(dtransform.dm) %*% V %*% dtransform.dm
    }

    if(!is.null(contrast)){
      referencepos <- match(reference, object$x)
      if(contrast == "difference"){
        dcontrast.dtransform <- diag(nX)
        dcontrast.dtransform[referencepos, ] <- -1
        dcontrast.dtransform[referencepos, referencepos] <- 0
        est <- est-est[referencepos]
      }
      if(contrast == "ratio"){
        dcontrast.dtransform <- diag(1 / est[referencepos], nrow=nX, ncol=nX)
        dcontrast.dtransform[referencepos, ] <- -est / est[referencepos]^2
        dcontrast.dtransform[referencepos, referencepos] <- 1
        est <- est / est[referencepos]
      }
      V <- t(dcontrast.dtransform) %*% V %*% dcontrast.dtransform
      V[referencepos, ] <- 0
      V[, referencepos] <- 0
    }

    se <-  sqrt(diag(V))
    
    if(CI.type == "plain"){
      theta <- qqq * se
      lower <- est - theta
      upper <- est + theta
    }
    if(CI.type == "log"){
      theta <- qqq * se / est
      lower <- est * exp(-theta)
      upper <- est * exp(theta)
    }
    temp <- as.matrix(cbind(est, se, lower, upper), nrow=length(est), ncol=4)
    dimnames(temp) <- list(object$x,
      c("Estimate", "Std. Error", paste("lower",CI.level), paste("upper",CI.level)))
    est.table[[j]] <- temp

  }
  if(is.factor(reference))
    reference <- as.character(reference)
  out <- c(object,
    list(est.table=est.table, tsum=t, transform=transform, contrast=contrast, reference=reference))
  class(out) <- "summary.stdCoxph"
  return(out)

}

print.summary.stdCoxph <- function(x, ...){

  nt <- length(x$tsum)
  for(j in 1:nt){

    cat("\nFormula: ")
    print(x$fit$formula)
    cat("Exposure: ", x$X, "\n")

    if(!is.null(x$transform))
      cat("Transform: ", x$transform,  "\n")
    if(!is.null(x$contrast)){
      cat("Reference level: ", x$X, "=", x$reference,  "\n")
      cat("Contrast: ", x$contrast,  "\n")
    }
    cat("Survival functions evaluated at t =", x$tsum[j], "\n")
    cat("\n")
    print(x$est.table[[j]], digits=3)
    cat("\n")

  }

}

plot.stdCoxph <- function(x, plot.CI=TRUE, CI.type="plain", CI.level=95,
  transform=NULL, contrast=NULL, reference=NULL, ...){

  object <- x
  x <- object$x

  dots <- list(...)

  xlab <- "t"

  if(is.factor(reference))
    reference <- as.character(reference)

  if(is.null(contrast)){
    if(is.null(transform))
      ylab <- expression(S(t))
    else{
      if(transform == "log")
        ylab <- expression(paste(log, "{", S(t), "}", sep=""))
      if(transform == "logit")
        ylab <- expression(paste(logit, "{", S(t), "}", sep=""))
      if(transform == "odds")
        ylab <- expression(paste(S(t), " / {", 1 - S(t), "}", sep=""))
    }
  }
  else{
    if(contrast == "difference"){
      if(is.null(transform))
        ylab <- c(bquote(paste(S(t), " - ", S[.(reference)](t))), expression())
      else{
        if(transform == "log")
          ylab <- c(bquote(paste(log, "{", S(t), "} - ", log, "{",
            S[.(reference)](t), "}", sep="")), expression())
        if(transform == "logit")
          ylab <- c(bquote(paste(logit, "{", S(t), "} - ", logit,
            "{", S[.(reference)](t), "}", sep="")), expression())
        if(transform == "odds")
          ylab <- c(bquote(paste(S(t), " / {", 1 - S(t), "} - ",
            S[.(reference)](t), " / {", 1 - S[.(reference)](t),
            "}", sep="")), expression())
      }
    }
    if(contrast == "ratio"){
      if(is.null(transform))
        ylab <- c(bquote(paste(S(t), "  /  ", S[.(reference)](t), sep="")), expression())
      else{
        if(transform == "log")
          ylab <- c(bquote(paste(log, "{", S(t), "}  /  ", log,
            "{", S[.(reference)](t), "}", sep="")), expression())
        if(transform == "logit")
          ylab <- c(bquote(paste(logit, "{", S(t), "}  /  ", logit,
            "{", S[.(reference)](t), "}", sep="")), expression())
        if(transform == "odds")
          ylab <- c(bquote(paste("[", S(t), " / {", 1 - S(t), "}]  /  [",
            S[.(reference)](t), " / {", 1 - S[.(reference)](t),
            "}]", sep="")), expression())
      }
    }
  }

  t <- object$t
  nt <- length(t)
  nX <- length(x)

  sum.obj <- summary(object=object, CI.type=CI.type, CI.level=CI.level,
    transform=transform, contrast=contrast, reference=reference)
  temp <- Reduce(f=rbind, x=sum.obj$est.table)
  est <- matrix(temp[, 1], nrow=nt, ncol=nX, byrow=TRUE)
  lower <- matrix(temp[, 3], nrow=nt, ncol=nX, byrow=TRUE)
  upper <- matrix(temp[, 4], nrow=nt, ncol=nX, byrow=TRUE)

  if(plot.CI)
    ylim <- c(min(lower), max(upper))
  else
    ylim <- c(min(est), max(est))

  args <- list(x=object$t, y=rep(0, length(t)), xlab=xlab, ylab=ylab,
    ylim=ylim, type="n")
  args[names(dots)] <- dots
  do.call("plot", args=args)

  legend <- NULL
  for(i in 1:nX){
    lines(t, est[, i], col=i)
    if(plot.CI){
      lines(t, upper[, i], col=i, lty="dashed")
      lines(t, lower[, i], col=i, lty="dashed")
    }
    temp <- as.character(x[i]) 
    legend <- c(legend, paste(object$X, "=", object$x[i])) 
  }
  if(is.na(match("ylim",names(args))))
    yl <- ylim[2]
  else
    yl <- args$ylim[2]
  legend(x=max(object$t), y=yl, legend=legend, lty=rep(1, length(x)),
      col=1:length(x), xjust=1)


}




