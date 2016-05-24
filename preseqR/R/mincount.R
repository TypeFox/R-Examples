## Two roots are the same if the difference is less than the PRECISION
PRECISION = 1e-3


### interpolating for species accumulation curve with minimum count
### ss step size
### n two-column histogram
### r minimum count
preseqR.interpolate.mincount <- function(ss, n, r=1)
{
  checking.hist(n)

  ## total individuals captured
  total.sample <- as.double(n[, 1] %*% n[, 2])
  N <- total.sample

  ## total species
  initial.distinct <- sum(as.numeric(n[, 2]))
  step.size <- as.double(ss)

  ## l is the number of sampled points for interpolation
  l <- N / step.size

  ## if the sample size is larger than the size of experiment or 
  ## the step size is too small, return NULL
  if (l < 1 || ss < 1 || r < 1)
    return()
  ## if the sample size is the size of the experiment
  ## count the number of species observed r or more times
  else if (l == 1) {
    index = which(n[, 1] >= r)
    result <- matrix(c(step.size, sum(as.numeric(n[index, 2]))), ncol = 2, byrow = FALSE)
    colnames(result) <- c('sample.size', 'interpolation')
    return(result)
  }

  ## explicitly calculating the expected species observed at least r times
  ## based on sampling without replacement
  ## see K.L Heck 1975
  ## N total individuals
  ## S the number of species
  ## size the size of the subsample
  expect.distinct <- function(n, N, size, S, r) {
    denom <- lchoose(N, size)
    p <- sapply(n[, 1], function(x) {
           sum(exp(lchoose(N - x, size - 0:(r-1)) + lchoose(x, 0:(r-1)) - denom))})
    return(S - p %*% n[, 2])
  }

  ## sample sizes
  x <- step.size * ( 1:l )

  ## calculate the number of distinct reads based on each sample size
  yield.estimates <- sapply(x, function(x) {
      expect.distinct(n, N, x, initial.distinct, r)})

  ## put size and yield together into a matrix
  result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}


### convert a continued fraction to a rational function
### the form of continued fractions e.g. a_0(t - 1) / (1 + a_1(t - 1))
### see the Supplementary of Daley, T., & Smith, A. D.(2013)
CF2RFA <- function(CF)
{
  poly.numer <- polynomial(1)
  poly.denom <- polynomial(1)
  cf <- CF$cf.coeffs
  for (i in rev(cf)) {
    tmp <- poly.numer
    poly.numer <- poly.numer + polynomial(c(0, i)) * poly.denom
    poly.denom <- tmp
  }
  ## according to the representation of the continued fraction
  ## the first item is a_0x/... not 1 + a_0x/...
  ## substract one from the result
  poly.numer <- poly.numer - poly.denom
  polylist(poly.numer, poly.denom)
}


### power series based on count frequencies starting from frequency j
### when j = 1, it is the power series expansion of E(S_1(t)) / t at t = 1
generating.ps <- function(n, j) {
  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]
  
  if (j >= length(hist.count)) {
    return(NULL)
  }

  ## shift to required count frequencies
  hist.count <- hist.count[j: length(hist.count)]

  PS.coeffs <- sum(n[ j:dim(n)[1], 2])
  change.sign = 0

  for (i in hist.count) {
    PS.coeffs <- c(PS.coeffs, (-1)^change.sign * i - PS.coeffs[length(PS.coeffs)])
    change.sign = change.sign + 1
  }

  ## truncate at coefficients where it is zero
  zero.index = which(PS.coeffs == 0)
  if (length(zero.index) > 0) {
    PS.coeffs[1:(min(zero.index) - 1)]
  } else {
    PS.coeffs
  }
}


### species accum curves based on parital fraction expansion of
### rational function approximation to E(S_1(t)) / t
### CHAO: the main function
preseqR.pf.mincount <- function(n, mt = 100, ss = NULL, 
                                max.extrapolation = NULL, r=1)
{
  # check the input format of the histogram
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## total individuals
  total.sample <- as.double(n[, 1] %*% n[, 2])

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- total.sample
    step.size <- ss
  } else if (ss < 1) {
    write("The step size is too small", stderr())
    return()
  } else {
    step.size <- ss
  }

  ## no interpolation if step.size is larger than the size of the iniital experiment
  ## set the starting sample size as the step.size
  if (step.size > total.sample) {
    yield.estimates <- NULL

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolating when sample size is no more than total sample size
      ## and setting the size of the sample for an initial extrapolation
      yield.estimates <- lapply(r, function(x) {preseqR.interpolate.mincount(step.size, n, x)})

      ## starting sample size for extrapolation
      starting.size <- ( floor(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## constructing the power series
  PS.coeffs <- generating.ps(n, 1)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## constrain the continued fraction approximation with even terms
  ## asymptotically ~ C / t
  mt <- min(mt, length(PS.coeffs))
  mt <- mt - (mt %% 2)
  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## check whether sample size is sufficient
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  ## construct a continued fraction approximation including as many as possible
  ## terms
  valid = FALSE
  DE = seq(mt, MIN_REQUIRED_TERMS, by=-2)
  for (de in DE) {
    ## continued fraction approximation to a power series
    out <- .C('c_PS2CF', as.integer(di), 
              as.integer(de), as.double(PS.coeffs[1:de]), 
              as.integer(length(PS.coeffs[1:de])),
              ps.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
              ps.coeffs.l = as.integer(0),
              cf.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
              cf.coeffs.l = as.integer(0),
              offset.coeffs =as.double(vector(mode='numeric',length=MAXLENGTH)),
              diagonal.idx = as.integer(0),
              degree = as.integer(0), is.valid = as.integer(0));

    if (out$is.valid) {
      ## pass results into R variables
      length(out$ps.coeffs) <- out$ps.coeffs.l
      length(out$cf.coeffs) <- out$cf.coeffs.l
      length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
      CF <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
                 out$diagonal.idx, out$degree)
      names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                     'degree')
      class(CF) <- 'CF'

      ## convert the continued fraction to the RFA and be compatible with 
      ## the constructed continued fraction in C extension
      RF <- CF2RFA(CF)
      RF[[1]] <- RF[[1]] / polynomial(c(0, 1))

      ## solving roots
      numer.roots <- solve(RF[[1]])
      denom.roots <- solve(RF[[2]])
      ## seperating roots by their real parts
      numer.roots.neg <- numer.roots[which(Re(numer.roots) < 0)]
      numer.roots.pos <- numer.roots[which(Re(numer.roots) >= 0)]
      denom.roots.neg <- denom.roots[which(Re(denom.roots) < 0)]
      denom.roots.pos <- denom.roots[which(Re(denom.roots) >= 0)]

      ## record roots in the numerator that are significantly similar to
      ## roots in the denominator
      tmp.roots <- c()

      ## simplify the rational function approximation
      ## two roots are same if the difference is less than the 
      ## predefined PRECISION
      if (length(numer.roots.pos) > 0) {
        for (i in 1:length(numer.roots.pos)) {
          if (length(denom.roots.pos) > 0) {
            d <- Mod(denom.roots.pos - numer.roots.pos[i])
            for (j in 1:length(d)) {
              if (d[j] < PRECISION) {
                denom.roots.pos <- denom.roots.pos[-j]
                tmp.roots <- c(tmp.roots, numer.roots.pos[i])
                break
              }
            }
          }
        }
      }

      ## roots in simplified RFA
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## convert roots from t - 1 to t
      roots <- denom.roots + 1
      ## max.range = max.extrapolation / total.sample
      ## pacman rule checking
      if (length(which(roots == 0)) || length(which(Re(roots) > 0))) {
        next
      } else {
        poly.numer <- as.function(poly.from.roots(numer.roots))
        l <- length(denom.roots)
        ## treat polynomials in the rational function to be monic
        ## the difference to the original RFA is a multiplier C

        ## c_i in the estimator
        coef <- sapply(1:l, function(x) {
          poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
        ## check whether the estimator is non-decreased
        deriv.f <- function(t) {
          Re(sapply(t, function(x) {-(coef*roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
        if (length(which( deriv.f(seq(0.05, as.double(max.extrapolation / total.sample), by=0.05)) < 0 ) != 0)) {
          next
        }
        ## the estimator passes the requirement
        valid = TRUE
        ## calculate the constant C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
        denom.roots <- denom.roots + 1
        coef <- coef * C
        mincount.f.elements <- list(coef, denom.roots)
        mincount.accum.curve.f <- lapply(r, function(x) {
            function(t) { Re(sapply(t, function(y) {coef %*% ( y / (y-denom.roots) )^x}))}})
        break
      }
    }
  }


  ## if the sample size is larger than max.extrapolation
  ## stop extrapolation
  ## MINOR.correction prevents machinary precision from biasing comparison
  ## result
  if (starting.size > (max.extrapolation + MINOR.correction))
  {
    if (length(yield.estimates) == 0) {
      return(NULL)
    } else {
      l <- 1:length(r)
      result <- lapply(l, function(x) {
          estimates <- yield.estimates[[x]]
          colnames(estimates) <- c("sample.size", paste("yield.estimates(r=", r[x], ")", sep=""))
          })
      return(result)
    }
  }

  ## valid is true if existing a RFA that satisfies pacman rule
  if (valid == FALSE) {
    return(NULL)
  }

  ## extrapolation for experiment with large sample size
  start <- starting.size / total.sample
  end <- (max.extrapolation + MINOR.correction) / total.sample
  step <- step.size / total.sample

  ## extrapolating for the general accumulation curves
  l = 1:length(r)
  extrap <- lapply(l, function(x) {
      mincount.accum.curve.f[[x]](seq(start, end, by=step))})

  result <- lapply(l, function(x) {
      ## combine results from interpolation and extrapolation
      estimates <- c(yield.estimates[[x]][, 2], extrap[[x]])
      index <- as.double(step.size) * (1: length(estimates))
      ## put index and estimated yields together into a two-colunm matrix
      estimates <- matrix(c(index, estimates), ncol = 2, byrow = FALSE)
      colnames(estimates) <- c("sample.size", paste("yield.estimates(r=", r[x], ")", sep=""))
      estimates
      })

  list(PF.elements=mincount.f.elements, yield.estimates=result)
}


## species accum curves based on parital fraction expansion
## the function is used whenever count frequency 1 is unavaible or the sample
## size is saturated
## using count frequencies starting from a given count frequency instead of 1
## when start.freq = 1, it is identical to the function preseqR.pf.mincount
## CHAO: save for a rainy day
general.preseqR.pf.mincount <- function(n, mt = 100, ss = NULL, 
                                        max.extrapolation = NULL, 
                                        r=1, start.freq=1)
{
  # check the input format of the histogram
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## total individuals
  total.sample <- as.double(n[, 1] %*% n[, 2])

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- total.sample
    step.size <- ss
  } else if (ss < 1) {
    write("The step size is too small", stderr())
    return()
  } else {
    step.size <- ss
  }

  ## no interpolation if step.size is larger than the size of the initial experiment
  if (step.size > total.sample) {
    yield.estimates <- vector(mode = 'numeric', length = 0)

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolation when sample size is no more than total sample size
      ## interpolate and set the size of the sample for an initial extrapolation
      out <- preseqR.interpolate.mincount(step.size, n, r)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( floor(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## power series for approximation 
  PS.coeffs <- generating.ps(n, start.freq)

  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }

  ## constrain the continued fraction approximation with even degree 
  ## asymptotically ~ C / t
  mt <- min(mt, length(PS.coeffs))
  mt <- mt - (mt %% 2)
  PS.coeffs <- PS.coeffs[ 1:mt ]

  ## check whether sample size is sufficient
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  ## indicator for existing an estimator satisfying the requirement
  valid = FALSE
  ## using as many terms as possible
  DE = seq(mt, MIN_REQUIRED_TERMS, by=-2)
  for (de in DE) {
    ## continued fraction approximation to a power series
    out <- .C('c_PS2CF', as.integer(di), 
              as.integer(de), as.double(PS.coeffs[1:de]), 
              as.integer(length(PS.coeffs[1:de])),
              ps.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
              ps.coeffs.l = as.integer(0),
              cf.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
              cf.coeffs.l = as.integer(0),
              offset.coeffs =as.double(vector(mode='numeric',length=MAXLENGTH)),
              diagonal.idx = as.integer(0),
              degree = as.integer(0), is.valid = as.integer(0));

    if (out$is.valid) {
      ## pass results into R variables
      length(out$ps.coeffs) <- out$ps.coeffs.l
      length(out$cf.coeffs) <- out$cf.coeffs.l
      length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
      CF <- list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, 
                 out$diagonal.idx, out$degree)
      names(CF) <- c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx',
                     'degree')
      class(CF) <- 'CF'

      ## convert the continued fraction to the RFA and be compatible with
      ## the constructed continued fraction in C extension
      RF <- CF2RFA(CF)
      RF[[1]] <- RF[[1]] / polynomial(c(0, 1))

      ## solving roots
      numer.roots <- solve(RF[[1]])
      denom.roots <- solve(RF[[2]])
      ## seperating roots by their real parts
      numer.roots.neg <- numer.roots[which(Re(numer.roots) < 0)]
      numer.roots.pos <- numer.roots[which(Re(numer.roots) >= 0)]
      denom.roots.neg <- denom.roots[which(Re(denom.roots) < 0)]
      denom.roots.pos <- denom.roots[which(Re(denom.roots) >= 0)]

      ## record roots in the numerator that are significantly similar to
      ## roots in the denominator
      tmp.roots <- c()

      ## simplify the rational function approximation
      ## two roots are considered to be identical if the difference is less than
      ## the predefined PRECISION
      if (length(numer.roots.pos) > 0) {
        for (i in 1:length(numer.roots.pos)) {
          if (length(denom.roots.pos) > 0) {
            d <- Mod(denom.roots.pos - numer.roots.pos[i])
            for (j in 1:length(d)) {
              if (d[j] < PRECISION) {
                denom.roots.pos <- denom.roots.pos[-j]
                tmp.roots <- c(tmp.roots, numer.roots.pos[i])
                break
              }
            }
          }
        }
      }

      ## roots in simplified RFA
      numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
      denom.roots <- c(denom.roots.neg, denom.roots.pos)

      ## convert roots from t - 1 to t for applying pacman rule
      roots <- denom.roots + 1
      ## assume 1 <= t <= max range
      ## max.range = max.extrapolation / total.sample
      ## pacman rule checking
      if (length(which(roots == 0)) || length(which(Re(roots) > 0 & 
          Mod(roots) / Re(roots) / 2 <= as.double(max.extrapolation / total.sample)))) {
        next
      } else {
        poly.numer <- as.function(poly.from.roots(numer.roots))
        l <- length(denom.roots)
        ## treat polynomials in the rational function to be monic
        ## the difference to the original RFA is a multiplier C

        ## c_i in the estimator
        coef <- sapply(1:l, function(x) {
          poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})

        ## modify the coefficients when using count frequencies not starting
        ## count frequency 1
        coef <- coef * (-denom.roots)^(start.freq - 1)
        ## check whether the estimator is non-decreased
        deriv.f <- function(t) {
          Re(sapply(t, function(x) {-(coef*roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
        if (length(which( deriv.f(seq(0.05, as.double(max.extrapolation / total.sample), by=0.05)) < 0 ) != 0)) {
          next
        }
        ## the estimator satisfies the requirement
        valid = TRUE
        ## the multiplier C
        C <- coef(RF[[1]])[length(coef(RF[[1]]))] / 
             coef(RF[[2]])[length(coef(RF[[2]]))]
        ## species accum curves with minimum count r
        ## using parital fraction expansion
        ## f is actually a function of "t -1" 
        f <- function(t) { 
          Re(sapply(t, function(x) {coef %*% ( (x+1)/(x-denom.roots) )^r})) * C}
        ## a function of "t"
        mincount.accum.curve.f <- function(x) {f(x - 1)}
        break
      }
    }
  }

  ## valid is true if existing a RFA that satisfies pacman rule
  if (valid == FALSE) {
    return(NULL)
  }

  ## if the sample size is larger than max.extrapolation
  ## stop extrapolation
  ## MINOR.correction prevents machinary precision from biasing comparison
  ## result
  if (starting.size > (max.extrapolation + MINOR.correction))
  {
    index <- as.double(step.size) * (1:length(yield.estimates))
    yield.estimates <- list(sample.size = index, yields = yield.estimates)
    result <- list(continued.fraction = CF, yield.estimates = yield.estimates)
    return(result)
  }

  ## extrapolation for experiment with large sample size
  start <- starting.size / total.sample
  end <- (max.extrapolation + MINOR.correction) / total.sample
  step <- step.size / total.sample

  ## extrapolating for the general accumulation curves
  extrap <- mincount.accum.curve.f(seq(start, end, by=step))

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, extrap)
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(estimator=mincount.accum.curve.f, estimates=yield.estimates)
  return(result)
}



### species accum curves based on parital fraction expansion of
### rational function approximation to E(S_1(t)) / t
### CHAO: the main function
preseqR.pf.mincount.bootstrap <- function(n, bootstrap.times = 20, mt = 100,
                                          ss = NULL, max.extrapolation = NULL, 
                                          conf = 0.95, r=1)
{
  # check the input format of the histogram
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## total individuals
  total.sample <- as.double(n[, 1] %*% n[, 2])

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- total.sample
    step.size <- ss
  } else if (ss < 1) {
    write("The step size is too small", stderr())
    return()
  } else {
    step.size <- ss
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }
  
  ## record second columns of resampled histograms                                  
  re.hist.second.col <- matrix(data = 0, nrow = length(n[, 1]),
  ncol = MULTINOMIAL.SAMPLE.TIMES) 

  ## the number of resampling times                                                 
  counter <- 0

  ## upperbound of times of iterations for bootstrapping
  upper.limit <- bootstrap.times/BOOTSTRAP.factor

  yield.estimates <- array(list(NULL), length(r))

  f <- function(x)
  {
    ## combine nonzero.index column and the second column to build a histogram
    ## table
    hist.table <- matrix(c(n[, 1], x), ncol = 2, byrow = FALSE)
    preseqR.pf.mincount(hist.table, mt=mt, ss=step.size, 
        max.extrapolation=max.extrapolation, r=r)
  }
  
  BOOTSTRAP.times = bootstrap.times

  while (bootstrap.times > 0) {
    ## do sampling with replacement
    ## re.hist.second.col saves the second columns of each resampled histogram
    re.hist.second.col <- rmultinom(MULTINOMIAL.SAMPLE.TIMES, sum(n[, 2]), n[, 2])

    ## estimate for each histogram
    out <- apply(re.hist.second.col, 2, f)

    ## eliminate NULL items in results
    out[sapply(out, is.null)] <- NULL
    ## extract yields estimation from each estimation result.
    yields <- lapply(out, function(x) x$yield.estimates)

    if ( length(yields) > 0 )
    {
      ## update sampling status
      success.times <- length(yields)
      bootstrap.times <- bootstrap.times - success.times
      for (i in 1:length(r)) {
        for (j in 1:length(yields)) {
          yield.estimates[[i]] <- cbind(yield.estimates[[i]], yields[[j]][[i]][, 2])
        }
      }
    }

    ## update sampling tmes
    counter <- counter + MULTINOMIAL.SAMPLE.TIMES
    if (counter > upper.limit)
      break;
  }

  ## enough successful sampling
  if (bootstrap.times <= 0) {
    
    if (BOOTSTRAP.times < 30) {
      warning("The confidence interval is not reliable because of insufficient iterations of bootstrapping.
  Set the variable bootstrap.times at least 30 in order to construct confidence intervals")}

    ## the number of sampled points for complexity curve
    n <- dim(yield.estimates[[1]])[1]

    ## sample sizes
    index <- as.double(step.size) * ( 1:n )
    result <- array(list(NULL), length(r))

    for (i in 1:length(r)) {

      # median values are used as complexity curve
      median.estimate <- apply(yield.estimates[[i]], 1, median)
      variance <- apply(yield.estimates[[i]], 1, var)

      # confidence interval based on lognormal
      if (conf <= 0 && conf >= 1)
        conf = 0.95
      C <- exp(qnorm((1 + conf) / 2.0) * sqrt(log(1.0 + variance / (median.estimate^2))))
      left.interval <- median.estimate/C
      right.interval <- median.estimate*C

      ## combine results and output a matrix
      result[[i]] <- matrix(c(index, median.estimate, left.interval, 
                            right.interval), ncol = 4, byrow = FALSE)
      lower.ci = sprintf('lower.%.2fCI', conf)
      upper.ci = sprintf('uppper.%.2fCI', conf)
      colnames(result[[i]]) <- c('sample.size', paste("yield.estimates(r=", r[i], ")", sep=""), lower.ci, upper.ci)
    }
    return(result)
  } else {
      write("fail to bootstrap!", stderr())
      return(NULL)
  }
}
