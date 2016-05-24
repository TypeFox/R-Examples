hdlm.fit <-
function(x, y, bootstrap=1, siglevel = 0.05, intercept=TRUE,
          alpha = 0.5, M = NULL, N = NULL, scale=TRUE, pval.method=c('mean', 'fdr', 'holm', 'QA'),
          ..., FUNCVFIT = NULL, FUNLM = NULL) {

  p <- ncol(x)
  n <- nrow(x)
  if(is.null(N)) N <- floor(n/2)
  if(is.null(M)) M <- floor((n - N) * 0.9)
  pval.method <- pval.method[[1]]

  # quietly load 'foreach':
  if (is.null(getDoParName())) {
    registerDoSEQ() 
  }

  NFOLDS <- 10
  if(n <= 200) NFOLDS <- max(c(3,floor(n/20)))

  if(is.null(FUNCVFIT)) FUNCVFIT <- function(x,y) {
    return(mod.cv.glmnet(x,y,alpha=alpha,standardize=scale,nfolds=NFOLDS))
  }
  if(is.null(FUNLM)) FUNLM <- lm

  hdlm.stat <- function(x, y, FUNCVFIT=FUNCVFIT, FUNLM=FUNLM, p=p, n=n, N=N) {
    INDEX <- sample(1:n,N)
    if(M == 0) INDEX <- 1:n
    #Estimate model:
    out <- FUNCVFIT(x[INDEX,],y[INDEX])
    if(M == 0) return(out)
    #lambda <- c(which(out$cvm <= min(out$cvm) + sd.off * out$cvsd[which.min(out$cvm)]), which.min(out$cvm))
    #lambda <- lambda[[1]]
    model <- which(out != 0)
 
    if(length(model) + intercept > M) {
      model <- order(abs(out), decreasing=TRUE)[1:M]
    }

    if(length(model) != 0 & intercept) {
      res <- rep(0,p+1)
      serr <- rep(0,p+1)
      out <- summary(FUNLM(y[-INDEX] ~ as.matrix(x[-INDEX,model])))
      if(sum(is.na(coef(out)[,1:2])) != 0 | nrow(coef(out)) != length(model)+1) {
          stop('FUNLM not reporting p-values; possible overfit model\nSee help pages for more information')
      }
      serr[c(1,model+1)] <- coef(out)[,2]
      res[c(1,model+1)]  <- coef(out)[,1]
    } else if(length(model) != 0 & !intercept)  {
      res <- rep(0,p)
      serr <- rep(0,p)
      out <- summary(FUNLM(y[-INDEX] ~ x[-INDEX,model] - 1))
      serr[c(model)] <- coef(out)[,2]
      res[c(model)] <- coef(out)[,1]
    } else if(length(model) == 0 & intercept)  {
      res <- rep(0,p+1)
      serr <- rep(0,p+1)
      out <- summary(FUNLM(y[-INDEX] ~ 1))
      serr[1] <- coef(out)[,2]
      res[1] <- coef(out)[,1]
    } else {
      res <- rep(0,p)
      serr <- rep(0,p)
    }
    index <- which(serr == 0)
    if(length(index) != 0 & length(index) != length(serr)) serr[index] <- min(serr[-index])
    return(cbind(res,serr))
  }

  # Calculate fdr pvalue given a vector of p-values 
  phi <- sum(1/(1:bootstrap))
  fdrpval <- function(v) {
    return(min(c(1,sort(v) * bootstrap * phi / (1:bootstrap))))
  }
  holmpval <- function(v) {
    return(min(c(1,min(v) * bootstrap)))
  }  
  QApval <- function(v) {
    gvals <- seq(0.05,1,length.out=floor(0.95*length(v)))
    v <- sort(v)[(length(v)-length(gvals)+1):length(v)]
    return(min(c(1,(1-log(0.05))*(v/gvals) )))
  }  

  #ht <- function(mu, u, v) {
  #  pvals <- (1-pnorm(abs(u - mu) / v))*2
  #  pvals[is.na(pvals) | v == 0] <- 1
  #  pvals[pvals > 1] <- 1
  #  return(fdrpval(pvals))
  #}

  #hs <- function(mu, u, v) {
  #  pvals <- (1-pnorm(abs(u - mu) / v))*2
  #  pvals[is.na(pvals) | v == 0] <- 1
  #  return(fdrpval(pvals))
  #}

  # If bootstrap = 1, make one run and collect results
  if(bootstrap == 1) {
    out <- hdlm.stat(x,y,FUNCVFIT=FUNCVFIT, FUNLM=FUNLM, p=p, n=n, N=N)
    if(M == 0) return(out)
    point_estimator <- out[,1]
    lb <- out[,1] - out[,2] * qnorm(1 - siglevel/2)
    ub <- out[,1] + out[,2] * qnorm(1 - siglevel/2)
    pvalue <- dnorm(point_estimator / out[,2])
    pvalue[is.na(pvalue)] <- 1
  } else {
    # If bootstrap != 1, calculate point estimates and SEs 'boostrap' times
    # Possible in parallel with 'foreach'

    output <- foreach(icount(bootstrap), .combine = "rbind", .inorder = FALSE) %dopar% {
        hdlm.stat(x,y,FUNCVFIT=FUNCVFIT, FUNLM=FUNLM, p=p, n=n, N=N)
    }
    #output <- matrix(0, nrow=bootstrap, ncol=p+intercept)
    vals <- matrix(output[,1], nrow=bootstrap, ncol=p+intercept, byrow=TRUE)
    ses <-  matrix(output[,2], nrow=bootstrap, ncol=p+intercept, byrow=TRUE)

    #vals <- matrix(0, nrow=bootstrap, ncol=p+intercept)
    #ses <- matrix(0,  nrow=bootstrap, ncol=p+intercept)
    #for(k in 1:bootstrap) {
    #  out <- output[[k]]
    #  vals[k,] <- out[,1]
    #  ses[k,] <- out[,2]
    #}

    # Calculate pvalues (matrix: variables by runs) from bootstrap estimates
    pvals <- (1-pnorm(abs(vals) / ses))*2
    pvals[is.na(pvals)] <- 1
    pvals[pvals > 1] <- 1

    # Collect pvalues across columns (runs) by use of fdrpval function above
    if(pval.method == 'median') {
      pvalue <- 2 * apply(pvals,2,median)
      if(sum(pvalue > 1) != 0) pvalue[pvalue > 1] <- 1
    } else if (pval.method == 'fdr') {
      pvalue <- apply(pvals,2,fdrpval)
    } else if (pval.method == 'holm') {
      pvalue <- apply(pvals,2,holmpval)
    } else {
      pvalue <- apply(pvals,2,QApval)
    }

    # Calculate fdr confidence intervals using function ht() from above
    lb <- rep(0, p + intercept)
    ub <- rep(0, p + intercept)
    point_estimator <- rep(0, p + intercept)
    for(j in 1:(p + intercept)) {
      if(sum(vals[,j] != 0) != 0) {
        #points <- seq(min(vals[,j]) - qnorm(1-siglevel/2) * max(abs(ses[,j])),
        #              max(vals[,j]) + qnorm(1-siglevel/2) * max(abs(ses[,j])),
        #              length.out=100)
        #out <- unlist(lapply(points, ht, vals[,j], ses[,j]))
        #interval <- points[which(out > siglevel)]
        #if(length(na.omit(interval)) == 0) interval <- range(points)
        #lb[j] <- min(interval)
        #ub[j] <- max(interval)
        #interval <- points[which(out == max(out))]
        #point_estimator[j] <- median(interval)
        lbounds <- sort(vals[,j] - qnorm(1-siglevel/(2*bootstrap)) * abs(ses[,j]))
        ubounds <- sort(vals[,j] + qnorm(1-siglevel/(2*bootstrap)) * abs(ses[,j]),decreasing=TRUE)
        temp <- floor((bootstrap + 1) / 2)
        lb[j] <- lbounds[temp]
        ub[j] <- ubounds[bootstrap - temp + 1]
        #lb[j] <- min(lbounds)
        #ub[j] <- max(ubounds)
        temp <- vals[,j]
        temp <- temp[temp >= lb[j] & temp <= ub[j]]
        point_estimator[j] <- median(temp, na.rm=TRUE)
        #point_estimator[j] <- (lb[j] + ub[j])/2 #median(vals[ lb[j] <= vals[,j] & ub[j] >= vals[j],j], na.rm=TRUE)
      }
    }


  }

  # Once we have a point estimator, calculated fitted values, residuals, and sigma hat
  # in the same way as one would for a non-penalized linear fit.
  if(intercept == TRUE) {
    fitted <- cbind(1,x) %*% point_estimator
  } else {
    fitted <- x %*% point_estimator    
  }
  resid <- fitted - y
  sigma_hat <- sd(as.numeric(resid))

  # Collect objects as a list and retrun to hdlm() function
  z <- list(coefficients=point_estimator, lower.bound=lb, upper.bound=ub, p.value=pvalue,
            effects=NULL, rank=c(n,p), fitted.values=fitted, assign=NULL, residuals=resid,
            sigma.hat=sigma_hat)

  return(z)

}

