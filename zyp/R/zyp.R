require(Kendall)

zyp.sen <- function(formula, dataframe) {
  zyp.slopediff <- function(i, xx, yy, n)
    (yy[1:(n-i)] - yy[(i+1):n]) / (xx[1:(n-i)] - xx[(i+1):n])

  if (missing(dataframe))
    dataframe <- environment(formula)
  term <- as.character(attr(terms(formula), "variables")[-1])
  x <- dataframe[[term[2]]]
  y <- dataframe[[term[1]]]
  n <- length(x)
  if (length(term) > 2) {
    stop("Only linear models are accepted")
  }
  
  slopes <- unlist(lapply(1:(n-1), zyp.slopediff, x, y, n))
  sni <- which(is.finite(slopes))
  slope <- median(slopes[sni])
  intercepts <- y - slope * x
  intercept <- median(intercepts)

  res <- list(coefficients=c(intercept, slope), slopes = slopes, intercepts = intercepts, rank=2, residuals=(y - slope * x + intercept), x=x, y=y)
  
  names(res$coefficients) = c("Intercept", term[2])
  class(res) = c("zyp", "lm")
  return(res)
}

confint.zyp <- function (object, parm, level = 0.95, ...) {
  res <- rep(0, 4)
  dim(res) <- c(2, 2)
  slopes <- sort(object$slopes)
  intercepts <- sort(object$intercepts)

  uquantile <- 1 - (1 - level) / 2
  zstat <- qnorm(uquantile)
  k <- Kendall(object$x, object$y)
  c.alpha <- sqrt(k$varS) * zstat
  n.prime <- length(slopes)
  isect.sd <- sqrt(var(intercepts))
  isect.mean <- mean(intercepts)
  
  idx <- c(round((n.prime - c.alpha) / 2), round((n.prime + c.alpha) / 2))
  
  rownames(res) <- names(object$coefficients)
  colnames(res) <- as.character(c((1 - level)/2, 1 - (1 - level)/2))
  res[2, ] <- slopes[idx]
  res[1, ] <- isect.mean + c(-isect.sd * zstat, isect.sd * zstat)
  return(res)
}

zyp.zhang <- function(y, x=1:length(y), conf.intervals=TRUE, preserve.range.for.sig.test=TRUE) {
  data <- as.numeric(as.vector(y))
  if(is.logical(x))
    stop("x cannot be of type 'logical' (perhaps you meant to specify conf.intervals?)")

  n <- length(data)
  t <- x

  ret <- c(lbound = NA, trend = NA, trendp = NA, ubound = NA,
           tau = NA, sig = NA, nruns = NA, autocor = NA, valid_frac = NA, linear = NA, intercept = NA)
  
  # Prewhiten the original series
  c <- acf(data,lag.max=1,plot=FALSE,na.action=na.pass)$acf[2]
  
  # Throw out data that is insufficient to compute the autocorrelation function
  if(is.na(c)) {
    return(ret)
  }

  if(c < 0.05) {
    y <- data
    yt <- t
  } else {
    y <- ifelse(rep(preserve.range.for.sig.test, n - 1), (data[2:n] - c * data[1:(n-1)]) / (1 - c), data[2:n] - c * data[1:(n-1)])
    yt <- t[1:(n-1)]
  }

  dmap <- which(!is.na(y))        
  ynm <- as.numeric(y[dmap])
  ytnm <- as.numeric(yt[dmap])
  
  # Throw out crappy data
  if(length(dmap) <= 3 | length(which(ynm != 0)) < 3 | length(dmap) / n < 0.1) {
    return(ret)
  }

  # Calculate Slope
  sen <- zyp.sen(ynm~ytnm)
  trend <- sen$coefficients[2]
        
  k <- 1
  c0 <- c
  if(c >= 0.05) {
    # Use repeated pre-whitening if there is autocorrelation
    trend0 <- trend
    while(k<500) {
      x <- data[1:n] - trend * t
      c <- acf(x,lag.max=1,plot=FALSE,na.action=na.pass)$acf[2]
      if(c < 0.05 && abs(c-c0) <= 0.0001) {
        break;
      }
      y <- (data[2:n] - c * data[1:(n-1)]) / (1 - c)
            
      dmap <- which(!is.na(y))
      ynm <- as.numeric(y[dmap])
      ytnm <- as.numeric(yt[dmap])
      sen <- zyp.sen(ynm~ytnm)
      trend <- sen$coefficients[2] # median(na.omit(sen$slopes))

      k <- k+1
            
      tpdiff <- abs((trend - trend0) / trend)
      if(!is.na(tpdiff) && tpdiff <= 0.001 && abs(c-c0)<=0.0001){
        break
      }
            
      trend0 <- trend
      c0 <- c
    }
  }
          
  Kend <- Kendall(ytnm,ynm)
  tau <- Kend[1]
  Bsig <- Kend[2]
  
  if(conf.intervals) {
    ci <- confint(sen)
  } else {
    ci <- matrix(rep(NA, 4), nrow=2, ncol=2)
  }
  
  ret <- c(lbound = as.numeric(ci[2, 1]), trend=as.numeric(trend), trendp=(as.numeric(trend) * n), ubound = as.numeric(ci[2, 2]),
           tau=as.numeric(tau), sig=as.numeric(Bsig), nruns=as.numeric(k), autocor=as.numeric(c0), valid_frac=as.numeric(length(dmap)/length(y)),
           linear=as.numeric(lm(data~t)$coefficients[2]), intercept=as.numeric(sen$coefficients[1]))
  return(ret)
}

zyp.yuepilon <- function(y, x=1:length(y), conf.intervals=TRUE, preserve.range.for.sig.test=TRUE) {
  dat <- as.numeric(as.vector(y))

  if(is.logical(x))
    stop("x cannot be of type 'logical' (perhaps you meant to specify conf.intervals?)")
  
  n <- length(dat)
  t <- x
  t.prime <- t[1:(n-1)]
  y <- dat

  ret <- c(lbound = NA, trend=NA, trendp=NA, ubound = NA,
           tau=NA, sig=NA, nruns=NA, autocor=NA, valid_frac=NA, linear=NA, intercept=NA)
  
 
  dmap <- which(!is.na(y))        
  ynm <- as.numeric(y[dmap])
  tnm <- as.numeric(t[dmap])

  # Throw out crappy data
  if(length(dmap) <= 3 | length(which(ynm != 0)) < 3 | length(dmap) / n < 0.1) {
    return(ret)
  }

  # Calculate Slope
  sen <- zyp.sen(ynm~tnm)
  trend <- sen$coefficients[2]

  # FIXME: ADD CHECK HERE FOR SMALL TREND
  
  # Remove trend
  xt.prime <- dat[1:n] - trend * t
  
  # Calculate AR(1)
  ac <- acf(xt.prime, lag.max=1, plot=FALSE, na.action=na.pass)$acf[2]

  # Throw out data that is insufficient to compute the autocorrelation function
  if(is.na(ac)) {
    return(ret)
  }

  yt.prime <- ifelse(rep(preserve.range.for.sig.test, n - 1), (xt.prime[2:n] - ac * xt.prime[1:(n-1)]) / (1 - ac), xt.prime[2:n] - ac * xt.prime[1:(n-1)])
  
  ## Add the trend back into the residual
  yt <- yt.prime[1:(n-1)] + trend * t.prime
  dmap.prime <- which(!is.na(yt))        
  ytnm <- as.numeric(yt[dmap.prime])
          
  # Calculate the Mann-Kendall test for significance on the blended series
  Kend <- Kendall(t.prime[dmap.prime],ytnm)
  tau <- Kend[1]
  Bsig <- Kend[2]

  if(conf.intervals) {
    ci <- confint(sen)
  } else {
    ci <- matrix(rep(NA, 4), nrow=2, ncol=2)
  }

  ret <- c(lbound = as.numeric(ci[2, 1]), trend=as.numeric(trend), trendp=as.numeric(trend) * n, ubound = as.numeric(ci[2, 2]),
           tau=as.numeric(tau), sig=as.numeric(Bsig), nruns=1, autocor=as.numeric(ac), valid_frac=as.numeric(length(dmap)/length(y)),
           linear=as.numeric(lm(dat~t)$coefficients[2]), intercept=as.numeric(sen$coefficients[1]))
           
  return(ret)
}

# Applies either a Yue Pilon or Zhang trend calculation to a vector of data
zyp.trend.vector <- function(y, x=1:length(y), method=c("yuepilon", "zhang"), conf.intervals=TRUE, preserve.range.for.sig.test=TRUE) {
  if(is.character(x))
    stop("x cannot be of type character (perhaps you meant to specify method?)")

  switch(match.arg(method),
         yuepilon = zyp.yuepilon(y, x, conf.intervals, preserve.range.for.sig.test),
         zhang    = zyp.zhang(y, x, conf.intervals, preserve.range.for.sig.test)
         )
}

# Applies either a Yue Pilon or Zhang trend calculation to a data frame where there is one station per line with no metadata (and each year/month is a column of data)
zyp.trend.dataframe <- function(indat, metadata.cols, method=c("yuepilon", "zhang"), conf.intervals=TRUE, preserve.range.for.sig.test=TRUE) {
  trend <- switch(match.arg(method),
           yuepilon = as.data.frame(t(apply(indat[, (1 + metadata.cols):ncol(indat)], 1, zyp.yuepilon, conf.intervals=conf.intervals, preserve.range.for.sig.test=preserve.range.for.sig.test))),
           zhang    = as.data.frame(t(apply(indat[, (1 + metadata.cols):ncol(indat)], 1, zyp.zhang, conf.intervals=conf.intervals, preserve.range.for.sig.test=preserve.range.for.sig.test))) )

  if(metadata.cols > 0) {
    trend <- cbind(indat[,1:metadata.cols], trend)
    # Preserve metadata column names
    names(trend)[1:metadata.cols] <- names(indat)[1:metadata.cols]
  }
  return (trend)
}

# Applies either a Yue Pilon or Zhang trend calculation to a file in the format "metadata_1,...,metadata_n,year1,year2,...", where n is the number of metadata columns and the timeseries extending across the row
zyp.trend.csv <- function(filename, output.filename, metadata.cols, method=c("yuepilon", "zhang"), conf.intervals=TRUE, csv.header=TRUE, preserve.range.for.sig.test=TRUE) {
  indat <- read.csv(filename, csv.header)
  write.csv(zyp.trend.dataframe(indat, metadata.cols, method, conf.intervals, preserve.range.for.sig.test), output.filename, row.names=FALSE)
}

    
