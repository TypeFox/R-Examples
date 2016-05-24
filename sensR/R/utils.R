rescale <-
  function(pc, pd, d.prime, std.err, 
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle")) 
{
  m <- match.call(expand.dots = FALSE)
  m[[1]] <- as.name("list")
  m <- eval.parent(m) # evaluate the *list* of arguments
  arg <- c("pc", "pd", "d.prime")
  isPresent <- sapply(arg, function(arg) !is.null(m[[arg]]))
  if(sum(isPresent) != 1)
    stop("One and only one of pc, pd and d.prime should be given")
  method <- match.arg(method)
  Pguess <- ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  par <- arg[isPresent]
  if(!is.null(se <- m$std.err)) {
    stopifnot(is.numeric(se) && length(se) == length(m[[par]]))
    stopifnot(all(se[!is.na(se)] > 0))
  }
  if(par == "pc") {
    pc <- m[[par]]
    stopifnot(is.numeric(pc) && all(pc >= 0) && all(pc <= 1))
    tooSmall <- pc < Pguess
    pc[tooSmall] <- Pguess
    pd <- pc2pd(pc, Pguess)
    d.prime <- psyinv(pc, method = method)
    if(!is.null(se)) {
      se.pc <- se
      se.pc[tooSmall] <- NA
      se.pd <- se.pc / (1 - Pguess)
      se.d.prime <- se.pc / psyderiv(d.prime, method = method)
    }
  }
  if(par == "pd") {
    pd <- m[[par]]
    stopifnot(is.numeric(pd) && all(pd >= 0) && all(pd <= 1))
    pc <- pd2pc(pd, Pguess)
    d.prime <- psyinv(pc, method = method)
    if(!is.null(se)) {
      se.pd <- se
      se.pc <- se.pd * (1 - Pguess)
      se.d.prime <- se.pc / psyderiv(d.prime, method = method)
    }
  }
  if(par == "d.prime") {
    stopifnot(is.numeric(d.prime) && all(d.prime >= 0))
    d.prime <- m[[par]]
    pc <- psyfun(d.prime, method = method)
    pd <- pc2pd(pc, Pguess)
    if(!is.null(se)) {
      se.d.prime <- se
      se.pc <- se * psyderiv(d.prime, method = method)
      se.pd <- se.pc / (1 - Pguess)
    } 
  }
  coef <- data.frame(pc = pc, pd = pd, d.prime = d.prime)
  res <- list(coefficients = coef)
  if(!is.null(se))
    res$std.err <- data.frame(pc = se.pc, pd = se.pd,
                              d.prime = se.d.prime)
  res$method <- method
  class(res) <- "rescale"
  return(res)
}

print.rescale <- function(x, digits = getOption("digits"), ...)
{
  cat(paste("\nEstimates for the", x$method, "protocol:\n", sep = " "))
  print(coef(x))
  if(!is.null(x$std.err)) {
    cat("\nStandard errors:\n")
    print(x$std.err)
  }
  return(invisible(x))
}

pc2pd <- function(pc, Pguess)
### Maps pc to pd

### arg: pc: numeric vector; 0 <= pc <= 1
###      Pguess: the guessing probability; numeric scalar,
###              0 <= pc <= 1
### res: pd: numeric vector; 0 <= pc <= 1
{
  stopifnot(is.numeric(Pguess) && length(Pguess) == 1 &&
            Pguess >= 0 && Pguess <= 1)
  stopifnot(is.numeric(pc) && all(pc >= 0) && all(pc <= 1))
  pd <- (pc - Pguess) / (1 - Pguess)
  pd[pc <= Pguess] <- 0
  names(pd) <- names(pc)
  return(pd)
}

pd2pc <- function(pd, Pguess) {
### Maps pd to pc
  
### arg: pd: numeric vector; 0 <= pc <= 1
###      Pguess: the guessing probability; numeric scalar,
###              0 <= pc <= 1
### res: pc: numeric vector; 0 <= pc <= 1
  stopifnot(is.numeric(Pguess) && length(Pguess) == 1 &&
            Pguess >= 0 && Pguess <= 1)
  stopifnot(is.numeric(pd) && all(pd >= 0) && all(pd <= 1))
  pc <- Pguess + pd * (1 - Pguess)
  names(pc) <- names(pd)
  return(pc)
}

psyfun <-
  function(d.prime,
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle")) 
### Maps d.prime to pc for sensory discrimination protocols
  
### arg: d.prime: non-negative numeric vector
### res: pc: numeric vector
{
  method <- match.arg(method)
  stopifnot(all(is.numeric(d.prime)) && all(d.prime >= 0))
  psyFun <- switch(method,
                   duotrio = duotrio()$linkinv,
                   tetrad = tetrad()$linkinv,
                   triangle = triangle()$linkinv,
                   twoAFC = twoAFC()$linkinv,
                   threeAFC = threeAFC()$linkinv)
  pc <- numeric(length(d.prime))
### Extreme cases are not handled well in the links, so we need: 
  OK <- d.prime < Inf
  if(sum(OK) > 0)
    pc[OK] <- psyFun(d.prime[OK])
  pc[!OK] <- 1
  names(pc) <- names(d.prime)
  return(pc)
}

psyinv <- function(pc, 
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle")) 
### Maps pc to d.prime for sensory discrimination protocols

### arg: pc: numeric vector; 0 <= pc <= 1
### res: d.prime: numeric vector
{
  method <- match.arg(method)
  stopifnot(all(is.numeric(pc)) && all(pc >= 0) && all(pc <= 1))
  psyInv <- switch(method,
                   duotrio = duotrio()$linkfun,
                   tetrad = tetrad()$linkfun,
                   triangle = triangle()$linkfun,
                   twoAFC = twoAFC()$linkfun,
                   threeAFC = threeAFC()$linkfun)
  d.prime <- numeric(length(pc))
### Extreme cases are not handled well in the links, so we need: 
  OK <- pc < 1
  if(sum(OK) > 0)
    d.prime[OK] <- psyInv(pc[OK])
  d.prime[!OK] <- Inf
  names(d.prime) <- names(pc)
  return(d.prime)
}

psyderiv <-
  function(d.prime, 
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle")) 
### Computes the derivative of the psychometric functions at some
### d.prime for sensory discrimination protocols.
  
### arg: d.prime: non-negative numeric vector
### res: pc: numeric vector
{
  method <- match.arg(method)
  stopifnot(all(is.numeric(d.prime)) && all(d.prime >= 0))
  psyDeriv <- switch(method,
                     duotrio = duotrio()$mu.eta,
                     tetrad = tetrad()$mu.eta,
                     triangle = triangle()$mu.eta,
                     twoAFC = twoAFC()$mu.eta,
                     threeAFC = threeAFC()$mu.eta)
  Deriv <- numeric(length(d.prime))
### Extreme cases are not handled well in the links, so we need: 
  OK <- d.prime > 0 && d.prime < Inf
  if(sum(OK) > 0)
    Deriv[OK] <- psyDeriv(d.prime[OK])
  Deriv[d.prime == 0] <- NA
  Deriv[d.prime == Inf] <- 0
  names(Deriv) <- names(d.prime)
  return(Deriv)
}

## findcr <-
##   function (sample.size, alpha = .05, p0 = .5, pd0 = 0,
##             type = c("difference", "similarity"))
## {
##   ## Find the critical value of a one-tailed binomial test.
##   type <- match.arg(type)
##   ss <- sample.size
##   if(ss != trunc(ss) | ss <= 0)
##     stop("'sample.size' has to be a positive integer")
##   if(alpha <= 0 | alpha >= 1)
##     stop("'alpha' has to be between zero and one")
##   if(p0 <= 0 | p0 >= 1)
##     stop("'p0' has to be between zero and one")
##   if(pd0 < 0 | pd0 > 1)
##     stop("'pd0' has to be between zero and one")
##   ## Core function:
##   i <- 0
##   if(type == "difference") {
##     while (1 - pbinom(i, ss, pd0 + p0*(1-pd0)) > alpha) i <- i + 1
##     i + 1
##   }
##   else {
##     while(pbinom(i, ss, pd0 + p0*(1-pd0)) < alpha) i <- i + 1
##     i - 1
##   }
## }

test.crit <-
  function(xcr, sample.size, p.correct = 0.5, alpha = 0.05, test)
### Is xcr the critical value of a one-tailed binomial test?
### Result: boolean

### OBS: there is deliberately no requirement that xcr should be
### positive or less than sample.size.
{  
  if(test %in% c("difference", "greater")) ## alternative is "greater"
    ((1 - pbinom(q = xcr - 1, size = sample.size, prob = p.correct) <= alpha) &&
     (1 - pbinom(q = xcr - 2, size = sample.size, prob = p.correct) > alpha))
  else if(test %in% c("similarity", "less")) ## alternative is "less"
    ((pbinom(q = xcr, size = sample.size, prob = p.correct) <= alpha) &&
     (pbinom(q = xcr + 1, size = sample.size, prob = p.correct) > alpha))
  else
    stop("unknown 'test' argument")
}

findcr <-
  function(sample.size, alpha = 0.05, p0 = 0.5, pd0 = 0,
           test = c("difference", "similarity"))
### Find the critical value of a one-tailed binomial
### test. "difference" means a "greater" alternative hypothesis and
### "similarity" means a "less" alternative hypothesis.

### FIXME: What should this function do/return if the critical value
### is larger than the sample size? Or when it is negative as can
### happen with similarity? Examples:
### (xcr <- findcr(sample.size = 1, p0 = psyfun(1, "twoAFC"))) ## 2
### (xcr <- findcr(sample.size = 1, test = "similarity")) ## -1
### This means that there is no number large/small enough for this
### sample size that could give a significant p-value. Maybe this
### should just be a deliberate feature.
{
  ## match and test arguments:
  test <- match.arg(test)
  ss <- sample.size
### FIXME: Does this test work as intented?
  if(ss != trunc(ss) | ss <= 0)
    stop("'sample.size' has to be a positive integer")
  if(alpha <= 0 | alpha >= 1)
    stop("'alpha' has to be between zero and one")
  if(p0 <= 0 | p0 >= 1)
    stop("'p0' has to be between zero and one")
  if(pd0 < 0 | pd0 > 1)
    stop("'pd0' has to be between zero and one")
  ## core function:
  pc <- pd2pc(pd0, p0)
  if(test == "difference") {
    crdiff <-  function(cr)
      1 - pbinom(q = cr - 1, size = ss, prob = pc) - alpha
    interval <- c(0, ss + 2) ## deliberately outside allowed range
  }
  else if(test == "similarity") {
    crdiff <- function(cr)
      pbinom(q = cr + 1, size = ss, prob = pc) - alpha
    interval <- c(-2, ss) ## deliberately outside allowed range
  }
  else ## should never occur
    stop("'test' not recognized") 
  xcr <- round(uniroot(crdiff, interval = interval)$root)
  ## is xcr the critical value?:
  is.crit <- test.crit(xcr = xcr, sample.size = ss, p.correct = pc,
                       alpha = alpha, test = test)
  if(is.crit) return(xcr)
  ## if uniroot fails, then do a simple search around the vicinity of
  ## the result from uniroot:
  max.iter <- 20 ## avoid infinite loop
  xcr <- delimit(xcr - 10, lower = -1)
  i <- 0
  if(test == "difference") {
    while(1 - pbinom(q = xcr + i, size = ss, prob = pc) > alpha) {
      if(i > max.iter || xcr + i > ss) break 
      i <- i + 1
    }
    xcr <- xcr + i + 1
  }
  if(test == "similarity") {
    while(pbinom(q = xcr + i, size = ss, prob = pc) < alpha) {
      if(i > max.iter || xcr + i > ss) break
      i <- i + 1
    }
    xcr <- xcr + i - 1
  }
  ## is xcr now the critical value?:
  is.crit <- test.crit(xcr = xcr, sample.size = ss, p.correct = pc,
                       alpha = alpha, test = test)
  if(is.crit) return(xcr)
  else stop("Failed to find critical value")
}

delimit <- function(x, lower, upper, set.na = FALSE)
### Sets the values of x < lower to lower and values of x > upper to
### upper. If set.na is TRUE values are set to NA. If both lower and
### upper are supplied, the (lower < upper) has to be TRUE.
{
  m <- match.call()
  m[[1]] <- as.name("list")
  m <- eval.parent(m)
  if(!is.null(m$lower) && !is.null(m$upper))
    stopifnot(m$lower < m$upper)
  if(!is.null(m$lower))
    x[x < m$lower] <- ifelse(set.na, NA, m$lower)
  if(!is.null(m$upper))
    x[x > m$upper] <- ifelse(set.na, NA, m$upper)
  return(x)
}

normalPvalue <-
### Computes the p-value for a statistic that follows a standard
### normal distribution under the null hypothesis.

### Arguments:
### statistic - a numerical (vector?)
### alternative - the type of alternative hypothesis

### Value:
### the p-value, possibly a vector.
  function(statistic, alternative = c("two.sided", "less", "greater"))
{
  alternative <- match.arg(alternative)
  stopifnot(all(is.finite(statistic)))
  p.value <-
    switch(alternative,
           "greater" = pnorm(statistic, lower.tail = FALSE),
           "less" = pnorm(statistic, lower.tail = TRUE),
           "two.sided" = 2 * pnorm(abs(statistic), lower.tail = FALSE))
  return(p.value)
}


## Do not partially match arguments.
## If possible give functions explicitly named  arguments - preferably
##   with default values.
## Value readability over speed.
## Value accuracy over speed.
## Use small functions with conceptual - easy-to-understand tasks.
## 
