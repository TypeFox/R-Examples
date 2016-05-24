BrierScoreDecomposition <- function(p, y, calibration=list(method="bin", bins=10), probs=NA, n.boot=0) {

  p <- as.vector(p)
  y <- as.vector(y)

  stopifnot(all(p >= 0), 
            all(p <= 1), 
            all(y==0 | y==1),
            "method" %in% names(calibration),
            length(y) == length(p))

  # estimate calibration function P(y=1|p)
  if (calibration[["method"]] == "bin") {
    # estimate P(y=1|p) by NumberOf(y=1 and p in bin i) / NumberOf(p in bin i)
    bins <- calibration[["bins"]]
    # define number of bins and bin breaks
    stopifnot(length(bins) > 0)
    if (length(bins) == 1) {
      nbins <- bins
      brx <- seq(0, 1, length.out = nbins + 1)
    } else {
      bins = sort(bins)
      stopifnot(bins[1] <= 0, bins[length(bins)] >= 1, 
                all(bins >= 0), all(bins <=1))
      nbins <- length(bins) - 1
      brx <- bins
    } 
    # estimate conditional event frequency
    p.hist <- hist(x=p, breaks=brx, plot=FALSE, 
                   include.lowest=TRUE)$counts
    cond.counts <- hist(x=p[y==1], breaks=brx, plot=FALSE,
                        include.lowest=TRUE)$counts
    cond.freq <- cond.counts / p.hist

    # match p and P(y=1|p)
    p.bin <- cut(x=p, breaks=brx, include.lowest=TRUE, labels=FALSE)
    cal <- cond.freq[ p.bin ]
  } else if (calibration[["method"]] == "logistic") {
    # estimate P(y=1|p) by a logistic regression model 
    w <- getOption("warn")
    options(warn=-1)
    cal <- glm(y~p, family="binomial")$fitted
    options(warn=w)
  } else {
    stop(paste("unknown calibration method:", calibration[["method"]]))
  }

  # estimate climatology
  clim <- mean(y, na.rm=TRUE)

  # calculate the components
  REL <- mean((p - cal)^2, na.rm=TRUE)
  RES <- mean((cal - clim)^2, na.rm=TRUE)
  UNC <- clim * (1 - clim)


  # quantiles of the sampling distribution 
  N <- length(y)
  boot.qntls <- NA
  if (!any(is.na(probs)) & n.boot > 1) {
    stopifnot(all(probs > 0 & probs < 1))
    probs <- sort(probs)
    boot <- t(replicate(n.boot, {

      inds <- sample(1:N, N, replace=TRUE)
      bp <- p[inds]
      by <- y[inds]
      
      # bootstrap calibration function 
      if (calibration[["method"]] == "bin") {
        bins <- calibration[["bins"]]
        stopifnot(length(bins) > 0)
        if (length(bins) == 1) {
          nbins <- bins
          brx <- seq(0, 1, length.out = nbins + 1)
        } else {
          bins = sort(bins)
          stopifnot(bins[1] <= 0, bins[length(bins)] >= 1, 
                    all(bins >= 0), all(bins <=1))
          nbins <- length(bins) - 1
          brx <- bins
        } 
        p.hist <- hist(x=bp, breaks=brx, plot=FALSE, 
                       include.lowest=TRUE)$counts
        cond.counts <- hist(x=bp[by==1], breaks=brx, plot=FALSE,
                            include.lowest=TRUE)$counts
        cond.freq <- cond.counts / p.hist
    
        # match p and P(y=1|p)
        p.bin <- cut(x=bp, breaks=brx, include.lowest=TRUE, labels=FALSE)
        cal <- cond.freq[ p.bin ]
      } else if (calibration[["method"]] == "logistic") {
        # estimate P(y=1|p) by a logistic regression model 
        w <- getOption("warn")
        options(warn=-1)
        cal <- glm(by~bp, family="binomial")$fitted
        options(warn=w)
      } 
      # estimate climatology
      clim <- mean(by, na.rm=TRUE)

      # calculate the components
      REL <- mean((bp - cal)^2, na.rm=TRUE)
      RES <- mean((cal - clim)^2, na.rm=TRUE)
      UNC <- clim * (1 - clim)
      c(REL=REL, RES=RES, UNC=UNC)
    }))
    boot.qntls <- apply(boot, 2, quantile, probs=probs)
    rownames(boot.qntls) <- paste(probs)
  
    ret <- list(c(REL=REL, RES=RES, UNC=UNC),BOOT.quantiles=boot.qntls)
  } else {
    ret <- c(REL=REL, RES=RES, UNC=UNC)
  }

  # return the decomposition
  return(ret)
}

