##Gets postpr output whatever method was used
getppout <- function(x) {
  if (x$method == "rejection") {
    l <- length(x$values)
    if(l > 0) {
      y <- table(x$values) / l
    } else {
      y <- NA
    }
  } else {
    y <- x$pred
  }
  return(y)
}

##Corrects for adjustment to prior by removing test data
priorcorrect <- function(pred, index, set) {
  out <- pred * table(index) / table(index[-set])
  out <- out / sum(out)
}

##Function to perform postpr given acceptance threshold rather than a tolerance
pp.eps <- function(eps, index, ...) {
  temp <- abc(tol=1, param=index, ...)$dist
  tol <- mean(temp <= eps)
  postpr(tol=tol, index=index, ...)
}

##Returns c(M[1,cols[1]], M[2,cols[2]], ...)
pull <- function(M, cols) { 
  out <- numeric(length(cols))
  for (i in seq(along.with=cols)) {
    out[i] <- M[i,cols[i]]
  }
  return(out)
}

##Returns quantiles of categorical distribution with probabilities pr
qcat <- function(u, pr) {
  pr <- as.numeric(pr)
  sapply(u, function(uu) {
    match(TRUE, uu <= cumsum(pr))
  })
}

mc.ci <-
function(raw, tol, eps, modname, modtrue, nbins=5, bintype=c("interval","quantile"), bw=FALSE, ...) {
  if (missing(tol) && missing(eps)) stop("Either tol or eps must be specified")
  usetol <- !missing(tol)
  modname <- as.character(modname)
  if (usetol) {
    x <- raw[raw$tol==tol,]
  } else {
    x <- raw[raw$eps==eps,]
  }
  mm <- which(names(x) == modname)
  if (nrow(x) == 0 | length(mm) == 0) return()
  inmod <- modtrue[x$testset] == modname
  temp <- matrix(nrow=nbins, ncol=4)
  colnames(temp) <- c("low", "high", "N", "x")
  if (bintype[1]=="interval") {
    temp[,1] <- c(-1,seq(1,nbins-1))/nbins ##First interval closed on left, all others open
    temp[,2] <- 1:nbins/nbins
  } else {
    temp2 <- quantile(x[,modname], probs=seq(0,1,length.out=nbins+1))
    temp[,1] <- temp2[1:nbins]
    temp[,2] <- temp2[-1]
  }
  for (i in 1:nbins) {
    zz <- which(x[,mm] > temp[i,1] & x[,mm] <= temp[i,2])
    temp[i,3] <- length(zz)
    temp[i,4] <- sum(inmod[zz])
  }
  if (bintype[1]=="inteval") temp[1,1] <- 0 ##-1 would mess up plot
  temp <- as.data.frame(temp)
  temp <- cbind(temp, alpha=1+temp$x, beta=1+temp$N-temp$x)
  if (bintype[1]=="interval") {
    plot.defaults <- list(xlim=c(0,1), ylim=c(0,1), xlab="Predicted", ylab="Observed", type='l')
  } else {
    plot.defaults <- list(xlim=c(temp2[1]*0.99,temp2[nbins+1]*1.01), ylim=c(0,1), xlab="Predicted", ylab="Observed", type='l')
  }
  ##Plot arguments defined in this way to allow user to override default choices through ... without errors
  plot.args <- modifyList(plot.defaults, list(...))
  plot.args <- modifyList(plot.args, list(x=c(0,1), y=c(0,1)))
  do.call(plot, plot.args)
  if (bw) {
    colours <- rep("black", nbins)
  } else {
    colours <- rep(c("red","blue"), ceiling(nbins/2))
  }
  for (i in 1:nbins) {
    lines(c(0.97*temp$low[i]+0.03*temp$high[i], 0.97*temp$high[i]+0.03*temp$low[i]), rep(temp$alpha[i]/(temp$alpha[i]+temp$beta[i]), 2), col=colours[i])
    lines(c(0.97*temp$low[i]+0.03*temp$high[i], 0.97*temp$high[i]+0.03*temp$low[i]), rep(qbeta(0.025, temp$alpha[i], temp$beta[i]), 2), col=colours[i], lty=3)
    lines(c(0.97*temp$low[i]+0.03*temp$high[i], 0.97*temp$high[i]+0.03*temp$low[i]), rep(qbeta(0.975, temp$alpha[i], temp$beta[i]), 2), col=colours[i], lty=3)
  }
}

covstats.mc <-
function(raw, index, diagnostics=c("freq", "loglik.binary", "loglik.multi"), nacc.min=20) {
  index.l <- levels(factor(index)) ##Same commands as postpr uses
  testsets <- unique(raw$testset)
  usetol <- "tol" %in% colnames(raw)
  tosplit <- ifelse(usetol, "tol", "eps")
  diag <- NULL
  ##COMPUTE DIAGNOSTICS
  if ("freq" %in% diagnostics) {
    ofreq <- table(index[testsets])
    myu <- matrix(runif(1E3*length(testsets)), nrow=1E3)
    diag <- dlply(raw, tosplit, ##I should be able to use ddply but it produces errors...
                  function(x) {
                    xx <- NULL
                    if (any(x$nacc<nacc.min) || sum(complete.cases(x)) < nrow(x)) {
                      xx <- rep(NA, ncol(x)-3)
                    }
                    x <- x[,-3] ##Remove nacc entries
                    if (is.null(xx)) {
                      xx <- numeric(ncol(x)-2)
                      for (i in 3:ncol(x)) {
                        pfreq <- apply(myu, 1, function(z){ sum(z<=x[,i]) })
                        xx[i-2] <- mean(ofreq[i-2] < pfreq) + 0.5*mean(ofreq[i-2] == pfreq) ##p-values for one-tailed tests
                      }
                    }
                    xx <- 1-2*abs(xx-0.5) ##p-values for two-tailed tests
                    return(data.frame(tol=x[1,2], parameter=colnames(x)[-(1:2)], pvalue=xx))
                  })
    diag <- rbind.fill(diag)
    diag <- cbind(diag, test="freq")
  }
  if ("loglik.binary" %in% diagnostics) {
    testmod <- index[testsets]
    myu <- matrix(runif(1E3*length(testsets)), nrow=1E3)
    temp <- dlply(raw, tosplit, ##I should be able to use ddply but it produces errors...
                  function(x) {
                    xx <- NULL
                    if (any(x$nacc<nacc.min) || sum(complete.cases(x)) < nrow(x)) {
                      xx <- rep(NA, ncol(x)-3)
                    }
                    x <- x[,-3] ##Remove nacc entries
                    if (is.null(xx)) {
                      xx <- numeric(ncol(x)-2)
                      for (i in 3:ncol(x)) {
                        y <- x[,i]
                        inmod <- (testmod == index.l[i-2])
                        oll <- sum(log(y[inmod])) + sum(log(1-y[!inmod]))
                        simmod <- apply(myu, 1, function(z){ z<=x[,i] })
                        pll <- apply(simmod, 2, function(z){
                          sum(log(y[z])) + sum(log(1-y[!z]))
                        })
                        xx[i-2] <- mean(oll < pll) + 0.5*mean(oll==pll) ##p-values for one-tailed tests
                      }
                    }
                    xx <- 1-2*abs(xx-0.5) ##p-values for two-tailed tests
                    return(data.frame(tol=x[1,2], parameter=colnames(x)[-(1:2)], pvalue=xx))
                  })
    temp <- rbind.fill(temp)
    temp <- cbind(temp, test="loglik.binary")
    if (is.null(diag)) {
      diag <- temp
    } else {
      diag <- rbind(diag, temp)
    }
  }
  if ("loglik.multi" %in% diagnostics) {
    testmod <- as.numeric(factor(index[testsets]))
    myu <- matrix(runif(1E3*length(testsets)), nrow=1E3)
    temp <- dlply(raw, tosplit, ##I should be able to use ddply but it produces errors...
                  function(x) {
                    xx <- NULL
                    if (any(x$nacc<nacc.min) || sum(complete.cases(x)) < nrow(x)) {
                      xx <- NA
                    }
                    x <- x[,-3] ##Remove nacc entries
                    if (is.null(xx)) {
                      oll <- sum(log(pull(x[,-(1:2)], testmod)))
                      simmod <- matrix(nrow=length(testsets), ncol=1E3)
                      for (i in 1:length(testsets)) {
                        simmod[i,] <- qcat(myu[,i], x[i,-(1:2)])
                      }
                      pll <- apply(simmod, 2, function(z){
                        sum(log(pull(x[,-(1:2)], z)))
                      })
                      xx <- mean(oll < pll) + 0.5*mean(oll==pll) ##p-value for one-tailed test
                    }
                    xx <- 1-2*abs(xx-0.5) ##p-value for two-tailed test
                    return(data.frame(tol=x[1,2], parameter=NA, pvalue=xx))
                  })
    temp <- rbind.fill(temp)
    temp <- cbind(temp, test="loglik.multi")
    if (is.null(diag)) {
      diag <- temp
    } else {
      diag <- rbind(diag, temp)
    }
  }
  if (!is.null(diag) && !usetol) {
    names(diag)[1] <- "eps"
  }
  return(diag)
}

cov.mc <-
function(index, sumstat, testsets, tol, eps, diagnostics=c(), multicore=FALSE, cores, method="rejection", nacc.min=20, ...) {
  ##CHECK INPUT
  if (!is.vector(index)) stop("index must be a vector")
  if (!is.data.frame(sumstat)) stop("sumstat must be a data frame")
  if (length(index) != nrow(sumstat)) stop("Number of rows in sumstat and length of index must be equal")
  if (missing(tol) && missing(eps)) stop("Either tol or eps must be specified")
  ##INITIALISE SOME VALUES
  n <- nrow(index)
  index.l <- levels(factor(index)) ##Same commands as postpr uses
  usetol <- !missing(tol)
  ##PREPARE INPUT VALUES
  if (usetol) {
    xx <- expand.grid(testset=testsets, tol=tol, stringsAsFactors=FALSE)
    tosplit <- .(tol)
  } else {
    xx <- expand.grid(testset=testsets, eps=eps, stringsAsFactors=FALSE)
    tosplit <- .(eps)
  }
  ##FUNCTION TO DO ONE ITERATION IE ABC FOR ONE INPUT CHOICE
  doit <- function(i){
    set <- xx$testset[i]
    if (usetol) {
      mytol <- xx$tol[i]
    } else {
      temp <- abc(target=sumstat[set,,drop=FALSE], param=index[-set], sumstat=sumstat[-set,,drop=FALSE], tol=1, method="rejection", ...)
      mytol <- mean(temp$dist <= xx$eps[i]) 
    }
    ##CHECK POST-PROCESSING POSSIBLE
    if (method != "rejection") {
      nacc <- floor(mytol * length(index))
      if (nacc < max(ncol(sumstat),length(index.l))) {
        out <- rep(NA, length(index.l))
        names(out) <- index.l
        out <- c(nacc=nacc, out)
        return(out)
      }
    }
    temp <- postpr(target=sumstat[set,,drop=FALSE], index=index[-set], sumstat=sumstat[-set,,drop=FALSE], tol=mytol, method=method, ...)
    pred <- getppout(temp)
    pred <- priorcorrect(pred, index, set)
    out <- c(nacc=length(temp$values), pred)
    return(out)
  }
  if (multicore) {
    raw <- mclapply(1:nrow(xx), doit, mc.cores=cores)
    raw <- simplify2array(raw)
  } else {
    raw <- sapply(1:nrow(xx), doit)
  }
  if (!is.matrix(raw)) {
    temp <- matrix(raw, nrow=1)
    rownames(temp) <- names(raw)[1]
    raw <- temp
  }
  raw <- cbind(xx,t(raw))
  ##RETURN OUTPUT
  return(list(raw=raw, diag=covstats.mc(raw, index, diagnostics, nacc.min)))
}
