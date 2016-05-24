### MAXLENGTH the allocate size for storing complexity curve
MAXLENGTH <- 10000000

### MULTINOMIAL.SAMPLE.TIMES number of random vectors to draw each time
MULTINOMIAL.SAMPLE.TIMES <- 11

### MINOR.correction a very small number to correct comparison result between
### two double type numbers when precisions can bias the result
MINOR.correction <- 1e-1

### BOOTSTRAP.factor the cut off ratio of success times / total bootstrap times
BOOTSTRAP.factor <- 0.4


### checking the input histogram in an appropariat format
checking.hist <- function(n)
{
  if (ncol(n)!=2 || is.numeric(n[,1])==FALSE || is.numeric(n[,2])==FALSE) {
    stop("Input must be a two-column matrix")
  }
  ## the first column is the frequencies of observed items
  freq <- n[, 1]

  ## the second column is the number of observed distinct items for each
  ## frequency
  number.items <- n[, 2]

  ## check whether frequencies are at least one and the histogram is sorted
  for (i in 1:length(freq))
    if (freq[i] <= 0 || freq[i] != floor(freq[i])) {
      stop("The first column must be positive integers!")
    } else if (number.items[i] < 0) {
      stop("The second column must be non negative")
    }
    else {
        if (i > 1 && freq[i - 1] >= freq[i])
          stop("The first column is not sorted in the ascending order")
    }

  return(n)
}


### calculate the value of the continued fraction approximation CF given the
### argument t
preseqR.rfa.estimate <- function(CF, t)
{
  if (class(CF) != "CF")
    return()
  if (t < 1) 
    stop("the t should be at least one!")
  
  ## the c function interface takes t - 1 as a parameter
  t <- t - 1
  ## call a c-encoded function c.calculate.continued.fraction
  out <- .C("c_calculate_continued_fraction",
            cf = as.double(CF$cf.coeffs),
            cf.l = as.integer(length(CF$cf.coeffs)),
            off = as.double(NULL),
            di = as.integer(0),
            de = as.integer(CF$degree),
            coordinate = as.double(t),
            result = as.double(0));

  ## return the calculated function value
  return(out$result * t)
}


### extrapolate given a histogram and a continued fraction
preseqR.extrapolate.distinct <- function(n, CF, start.size = NULL,
     step.size = NULL, max.size = NULL)
{
  checking.hist(n)
  ## check CF is a continued fraction with CF attribute
  if (class(CF) != "CF")
    return(NULL)

  ## parameters for calling the c-encode function c_extrapolate_distinct
  cf.coeffs <- as.double(CF$cf.coeffs)
  cf.coeffs.l <- as.integer(length(CF$cf.coeffs))
  offset.coeffs <- as.double(NULL)
  di <- as.integer(0)
  de <- as.integer(CF$degree)

  total.sample <- floor(n[, 1] %*% n[, 2])

  ## set start.size, step.size, max.size if they are not defined by user
  if (is.null(start.size))
    start.size <- 0
  if (is.null(max.size)) {
    ## 100 is a magic number
    max.size <- 100
  }
  if (start.size > max.size)
  {
    write("start position has already beyond the maximum prediction", stderr())
    return(NULL)
  }
  if (is.null(step.size))
    step.size <- 1

  ## allocate memory to store extrapolation results
  ## first "c.extrapolate.distinct" stores the observed number of distinct
  ## molecules into estimate, then it stores the extrapolation values
  ## thus the allocated memory size is 1 plus the size of extrapolation values,
  ## which is (max.size - start.size) / step.size) + 1
  extrap.size <- floor( (max.size - start.size) / step.size ) + 1

  out <- .C("c_extrapolate_distinct", cf.coeffs, cf.coeffs.l, offset.coeffs,
            di, de, as.double(start.size), 
            as.double(step.size), as.double(max.size), 
            estimate = as.double(vector(mode = 'numeric', extrap.size)),
            estimate.l = as.integer(0));

  initial_sum <- floor(sum(n[, 2]))
  extrapolation <- out$estimate[ 1:out$estimate.l ] + initial_sum

  ## sample size vector for extrapolation
  sample.size <- total.sample * (start.size + step.size*((1:length(extrapolation)) - 1)) + 
                 total.sample
  ## estimation should be conservative
  sample.size <- ceiling(sample.size)

  ## put sample.size and extrapolation results together into a matrix
  result <- matrix(c(sample.size, extrapolation), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'extrapolation')
  return(result)
}


## sampling without replacement
## n frequencies counts
nonreplace.sampling <- function(size, n)
{
  ## make sure frequencies are integers
  n[, 2] <- floor(n[, 2])
  ## the number of distinct items
  distinct <- sum(n[, 2])

  ## identifier for each distinct item
  ind <- 1:distinct

  ## the size of each read in the library
  N <- rep(n[, 1], n[, 2])

  ## construct a sample space X 
  ## the whole library represents by its indexes. If a read presents t
  ## times in the library, its indexes presents t times in X
  X <- rep(ind, N)

  return(sample(X, size, replace = FALSE))
}


## sampling without replacement
## input frequencies counts; output subsample as a frequencies counts
preseqR.nonreplace.sampling <- function(size, n)
{
  ## check the input histogram file
  checking.hist(n)
  ## sub sampling
  X <- nonreplace.sampling(size, n)
  ## record the freq of each sampled species
  freq.counts <- hist(X, breaks=0:max(X), plot=FALSE)$count
  ## frequencies counts; frequency 0 excluded
  T <- hist(freq.counts, breaks=-1:max(freq.counts), plot=FALSE)$counts[-1]
  matrix(c(which(T != 0), T[which(T != 0)]), byrow = FALSE, ncol=2)
}


### interpolate when the sample size is no more than the size of
### the initial experiment
preseqR.interpolate.distinct <- function(ss, n)
{
  checking.hist(n)

  ## calculate total number of sample
  total.sample <- n[, 1] %*% n[, 2]
  N <- floor(total.sample)

  initial.distinct <- sum(as.numeric(n[, 2]))
  ## the total individuals captured
  step.size <- as.double(ss)

  ## l is the number of interpolation points
  l <- as.integer(N / step.size)

  ## if the sample size is larger than the size of experiment or 
  ## the step size is too small, return NULL
  if (l == 0 || ss < 1)
    return()

  ## explicit calculating the expectation for sampling without replacement
  ## see K.L Heck 1975
  ## N is the size of population;   
  ## S is the number of unique species
  ## x is the size of the sub sample
  expect.distinct <- function(n, N, x, S) {
    denom <- lchoose(N, x)
    numer <- lchoose(N - n[, 1], x)
    p <- exp(numer - denom)
    return(S - p %*% n[, 2])
  }

  ## sample size vector
  x <- step.size * ( 1:l )

  ## calculate the number of distinct reads based on each sample size
  yield.estimates <- sapply(x, function(x) expect.distinct(n, N, x, initial.distinct))

  ## put sample.size and yield.estimates together into a matrix
  result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(result) <- c('sample.size', 'interpolation')

  return(result)
}

### check the goodness of the sample based on good Good & Toulmin's model
goodtoulmin.2x.extrap <- function(n)
{
  return((-1)^(n[, 1] + 1) %*% n[, 2])
}


### construct a rational function approximation given a frequencies of count
### data
### mt = max_terms, 
preseqR.rfa.curve <- function(n, mt = 100, ss = NULL,
                              max.extrapolation = NULL, asym.linear=FALSE)
{
  checking.hist(n)
  ## setting the diagonal value
  di = 0
  ## minimum required number of terms of power series in order to construct
  ## a continued fraction approximation
  MIN_REQUIRED_TERMS <- 4

  ## calculate total number of sample
  total.sample <- n[, 1] %*% n[, 2]
  total.sample <- floor(total.sample)

  ## set step.size as the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- floor(total.sample)
    step.size <- ss
  } else if (ss < 1) {
    write("step size is should be at least one", stderr())
    return(NULL)
  } else {
    step.size <- floor(ss)
  }

  ## no interpolation if step.size is larger than the size of experiment
  ## set the starting sample size as the step.size
  if (step.size > total.sample) {
    yield.estimates <- vector(mode = 'numeric', length = 0)

    ## starting sample size for extrapolation
    starting.size <- step.size
  } else {
      ## interpolation when sample size is no more than total sample size
      ## interpolate and set the size of sample for initial extrapolation
      out <- preseqR.interpolate.distinct(step.size, n)
      yield.estimates <- out[, 2]

      ## starting sample size for extrapolation
      starting.size <- ( as.integer(total.sample/step.size) + 1 )*step.size
  }

  if (is.null(max.extrapolation)) {
    ## extrapolation 100 times if it is undefined; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]
  ## only use non zeros items in histogram from begining up to the first zero
  counts.before.first.zero = 1
  while (as.integer(counts.before.first.zero) <= length(hist.count) &&
         hist.count[counts.before.first.zero] != 0)
    counts.before.first.zero <- counts.before.first.zero + 1

  ## constrain the continued fraction approximation with even degree 
  ## conservatively estimates
  mt <- min(mt, counts.before.first.zero - 1)
  if (asym.linear == TRUE && (mt %% 2 == 0)) {
    mt = mt - 1
  } else {
    mt = mt - (mt %% 2)
  }

  ## pre-check to make sure the sample is good for prediction
  if (mt < MIN_REQUIRED_TERMS)
  {
    m <- paste("max count before zero is les than min required count (4)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }
  if(goodtoulmin.2x.extrap(n) < 0.0)
  {
    m <- paste("Library expected to saturate in doubling of size",
               " unable to extrapolate", sep = ',')
    write(m, stderr())
    return(NULL)
  }

  ## adjust the format of count vector of the histogram in order to
  ## call the c-encoded function
  hist.count <- c(0, hist.count)

  ## allocate spaces to store constructed continued fraction approximation
  ## construct a continued fraction approximation with minimum degree
  out <- .C('c_continued_fraction_estimate', as.double(hist.count), 
            as.integer(length(hist.count)), as.integer(di), as.integer(mt),
            ps.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
            ps.coeffs.l = as.integer(0),
            cf.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
            cf.coeffs.l = as.integer(0),
            offset.coeffs =as.double(vector(mode='numeric',length=MAXLENGTH)),
            diagonal.idx = as.integer(0),
            degree = as.integer(0),
            is.valid = as.integer(0));

  if (!out$is.valid)
  {
    return(NULL)
  }

  ## restore the hist.count into the R coded format
  hist.count <- hist.count[-1]

  ## pass results into R variables
  length(out$ps.coeffs) <- out$ps.coeffs.l
  length(out$cf.coeffs) <- out$cf.coeffs.l
  length(out$offset.coeffs) <- as.integer(abs(out$diagonal.idx))
  CF <- list(out$ps.coeffs, out$cf.coeffs, out$degree)
  names(CF) <- c('ps.coeffs', 'cf.coeffs', 'degree')
  class(CF) <- 'CF'

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
  start <- ( starting.size - total.sample ) / total.sample
  end <- (max.extrapolation + MINOR.correction - total.sample) / total.sample
  step <- step.size / total.sample
  res <- preseqR.extrapolate.distinct(n, CF, start, step, end)

  ## combine results from interpolation/extrapolation
  yield.estimates <- c(yield.estimates, res[, 2])
  index <- as.double(step.size) * (1: length(yield.estimates))

  ## put index and estimated yields together into a two-colunm matrix
  yield.estimates <- matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE)
  colnames(yield.estimates) <- c('sample.size', 'yield.estimate')

  result <- list(continued.fraction = CF, estimates = yield.estimates)
  return(result)
}


### generate complexity curve through bootstrapping the histogram
preseqR.rfa.species.accum.curve <- function(
    n, bootstrap.times = 20, mt = 100, ss = NULL,
    max.extrapolation = NULL, conf = 0.95, asym.linear=FALSE)
{
  checking.hist(n)
  n[, 2] <- floor(n[, 2])
 
  ## setting the diagonal value
  di = 0

  ## calculate the total number of sample
  total.sample <- n[, 1] %*% n[, 2]

  ## set the step.size to the size of the initial experiment if it is undefined
  if (is.null(ss)) {
    ss <- total.sample
    step.size <- ss
  } else if (ss < 1) {
    write("step size should be at least one", stderr())
    return(NULL)
  } else {
    step.size <- floor(ss)
  }

  ## set the maximum extrapolation size if it is undefined
  if (is.null(max.extrapolation)) {

    ## extrapolation 100 times; 100 is a magic number
    max.extrapolation <- 100*total.sample
  }

  ## record second columns of resampled histograms
  re.hist.second.col <- matrix(data = 0, nrow = length(n[, 1]),
                               ncol = MULTINOMIAL.SAMPLE.TIMES)

  ## the number of resampling times
  counter <- 0

  yield.estimates <- vector(mode = "numeric", length = 0)

  ## upperbound of times of iterations for bootstrapping
  upper.limit <- bootstrap.times/BOOTSTRAP.factor

  f <- function(x)
  {
    ## combine nonzero.index column and the second column to build a histogram
    ## table
    hist.table <- matrix(c(n[, 1], x), ncol = 2, byrow = FALSE)
    preseqR.rfa.curve(hist.table, mt, step.size, max.extrapolation, asym.linear=asym.linear)
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
    yields <- sapply(out, function(x) x$estimates[, 2])

    if ( length(yields) > 0 )
    {
      ## update sampling status
      success.times <- dim(yields)[2]
      bootstrap.times <- bootstrap.times - success.times
      yield.estimates <- cbind(yield.estimates, yields)
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
    n <- dim(yield.estimates)[1]

    ## sample sizes
    index <- as.double(step.size) * ( 1:n )

    # median values are used as complexity curve
    median.estimate <- apply(yield.estimates, 1, median)
    variance <- apply(yield.estimates, 1, var)

    # confidence interval based on lognormal
    if (conf <= 0 && conf >= 1)
      conf = 0.95
    C <- exp(qnorm((1 + conf) / 2.0) * sqrt(log(1.0 + variance / (median.estimate^2))))
    left.interval <- median.estimate/C
    right.interval <- median.estimate*C

    ## combine results and output a matrix
    result <- matrix(c(index, median.estimate, left.interval, right.interval),
                    ncol = 4, byrow = FALSE)
    lower.ci = sprintf('lower.%.2fCI', conf)
    upper.ci = sprintf('uppper.%.2fCI', conf)
    colnames(result) <- c('sample.size', 'yield.estimate', lower.ci, upper.ci)
    return(result)
  } else {
      write("fail to bootstrap!", stderr())
      return(NULL)
  }
}

print.CF <- function(x, digit = 4, ...)
{
  s <- "CONTINUED FRACTION APPROXIMATION :\n\n"

  ## print the degree of the continued fraction approximation
  s <- paste(s, "DEGREE\t", toString(x$degree), "\n\n", sep = '')

  ## print the diagonal value
  ## s <- paste(s, "DIAGONAL VALUE\t", toString(x$diagonal.idx), "\n\n", sep = '')

  ## print the coefficients depending on the value of diagonal value
  s <- paste(s, "COEFFICIENTS:\n\n", sep = '')

  ##  di <- abs(x$diagonal.idx)
  di <- 0
    
  ## the function to print a coefficient
  ## S is the set of coefficients
  ## shift is the difference between the index of S and the index of 
  ## coefficients of a continued fraction approximation
  f <- function(index, S, shift)
  {
    s <- formatC(S[index], digit, format = 'f')
    s <- paste('a_', toString(index + shift), ' = ', s, sep = '')
  }
 
  ## print offset values if any
  if (di > 0)
  {
    index <- 1:di
    dim(index) <- di

    tmp <- apply( index, 1, function(t) f(t, x$offset.coeffs, -1) )
    tmp <- paste(tmp, collapse = '\n')
    s <- paste(s, tmp, '\n', sep = '')
  }

  ## print coeffients if any
  if (length(x$cf.coeffs) > 0)
  {
    index <- 1:length(x$cf.coeffs)
    dim(index) <- length(x$cf.coeffs)

    tmp <- apply( index, 1, function(t) f(t, x$cf.coeffs, di - 1) )
    tmp <- paste(tmp, collapse = '\n')
    s <- paste(s, tmp, '\n', sep = '')
  }
  cat(s)
  invisible(x)
}
