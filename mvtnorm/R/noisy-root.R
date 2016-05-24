## calculations of equi-coordinate quantiles for pmvt

wgts <- function(typred, ttarg){
  dists <- abs(typred-ttarg)
  if(sum(dists < 1 & dists > 0.0001) < 3){
    return(rep(1, length(typred)))
  } 
  (1-dists^3)^3*(dists < 1)
}

predict_with_se <- function(fit, newdata){
  pred <- try(predict(fit, newdata, se.fit=TRUE), silent=TRUE)
  if(inherits(pred, "try-error")){ ## problem with se calculation, use "large" se
    pred0 <- predict(fit, newdata)
    pred <- list(fit=pred0, se.fit=rep(0.5,length(pred0)))
  }
  pred
}

stop_crit <- function(ptol, fit, xest, ttarg, level=0.05){
  pred <- predict_with_se(fit, data.frame(x=xest))
  crit <- qt(1-level/2, df=fit$df.residual)
  cond1 <- pred$fit[1]+crit*pred$se.fit[1] < ttarg+ptol
  cond2 <- pred$fit[1]-crit*pred$se.fit[1] > ttarg-ptol
  return(list(stop=cond1 & cond2))
}

get_new_points <- function(fit, xest, targ, level=0.05){
  pred <- predict_with_se(fit, data.frame(x=xest))
  cf <- coef(fit)
  xLB <- (pred$fit-pred$se.fit-cf[1])/cf[2]
  xUB <- (pred$fit+pred$se.fit-cf[1])/cf[2]
  c(xLB, xUB)
}

get_est <- function(fit, ttarg){
  cf <- coef(fit)
  as.numeric(ttarg-cf[1])/cf[2]
}

fitlr <- function(x, y, weights, link){
  ind <- weights > 0
  weights <- weights[ind]
  x <- x[ind]  
  y <- y[ind]
  fit <- glm(y~x, weights = weights,
             family=quasi(link))
  return(fit)
}

sanitize_y <- function(y){
  eps <- .Machine$double.eps
  pmax(pmin(1-eps,y),eps)
}

get_quant_loclin <- function(func, targ, interval, 
                             ptol=0.005, maxiter = 500, 
                             link = c("probit", "cauchit"),
                             verbose = FALSE, ...){
  ## input argument checks
  if(interval[2] <= interval[1])
    stop("first entry of interval needs to be smaller than second")
  link <- match.arg(link)
  ## set up output vectors
  x <- y <- numeric(4+maxiter)
  xest <- numeric(1+maxiter)
  ## first three evaluations
  x[1:3] <- c(interval, mean(interval))
  res <- lapply(x[1:3], function(x) func(x, ...))
  y[1:3] <- sapply(res, function(x) sanitize_y(x))
  ## check for non-monotone function
  if(y[1] > y[2])
    stop("func does not appear to be monotone")
  if(link == "probit"){
    cdf <- pnorm
    quf <- qnorm
  } else {
    cdf <- pcauchy
    quf <- qcauchy
  }
  ttarg <- quf(targ) # transform target as well
  
  ## start iterating
  ycur <- y[1:3];xcur <- x[1:3]
  fit <- fitlr(xcur, ycur, weights=rep(1,3), link=link)
  xest[1] <- get_est(fit, ttarg)
  x[4] <- xest[1]
  y[4] <- sanitize_y(func(x[4], ...))
  count <- 4

  for(i in 1:maxiter){
    ind <- 1:count
    xcur <- x[ind];ycur <- y[ind]
    tycurpred <- predict(fit, newdata=data.frame(x=xcur))
    weights <- wgts(tycurpred, ttarg)
    fit <- fitlr(xcur, ycur, weights=weights, link=link)
    stp <- stop_crit(ptol, fit, xest[i], ttarg)
    if(verbose & i > 1){
      txt <- sprintf("Iteration %i, Estimate: %f, LB: %f (%f), UB: %f (%f)\n",
                     i, xest[i],
                     x[count-1], y[count-1],
                     x[count], y[count])
      cat(txt)
    }
    if(stp$stop)
      break
    xest[i+1] <- get_est(fit, ttarg)
    new_x <- get_new_points(fit, xest[i+1], targ)
    new_y <- lapply(new_x, function(x) func(x, ...))
    ind <- (count+1):(count+2)
    x[ind] <- new_x
    y[ind] <- sanitize_y(c(new_y[[1]], new_y[[2]]))
    count <- count+2
  }
  f.quant <- predict(fit, newdata=data.frame(x=xest[i]))
  out <- list(quantile = xest[i],
              f.quantile = cdf(as.numeric(f.quant)))
  if(i == maxiter){
    attr(out, "message") <- "Maximum number of iterations reached without sufficient accuracy"
  } else {
    attr(out, "message") <- "Normal Completion"
  }
  out
}

