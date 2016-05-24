##################################################################
### Power functions for the exact binomial power and the normal-based
### approximations thereto.
##################################################################
normalPwr <-
  function(pdA, pd0 = 0, sample.size, alpha = 0.05, pGuess = 1/2,
           test = c("difference", "similarity"), continuity = FALSE)
### Normal based approximation to the exact binomial power function
### optionally using a continuity correction.
{
  test <- match.arg(test)
  ss <- sample.size
  ## get
  pc0 <- pd2pc(pd=pd0, Pguess=pGuess)
  pcA <- pd2pc(pd=pdA, Pguess=pGuess)
  sigma0 <- sqrt(pc0*(1 - pc0)*ss)
  sigmaA <- sqrt(pcA*(1 - pcA)*ss)
  if(test == "difference") {
    lambda <- (qnorm(1 - alpha) * sigma0 + ss * (pc0 - pcA)) / sigmaA
    if(continuity) lambda <- lambda + .5/sigmaA
    pwr <- pnorm(lambda, lower.tail = FALSE)
  }
  else if(test == "similarity") {
    lambda <- (qnorm(alpha) * sigma0 + ss * (pc0 - pcA)) / sigmaA
    if(continuity) lambda <- lambda - .5/sigmaA
    pwr <- pnorm(lambda, lower.tail = TRUE)
  }
  else
    stop("'test' not recognized")
  return(as.vector(pwr))
}

pdPwr <-
  function(pdA, pcA, pd0 = 0, sample.size, alpha = 0.05, pGuess = 1/2,
           test = c("difference", "similarity"), crit=NULL)
### Simple, fast version of discrimPwr that does no testing of
### arguments and only provides power for the exact test.
{
  ## match and test arguments:
  test <- match.arg(test)
  ss <- sample.size
  ## Get pc from pdA and pGuess:
  if(missing(pcA))  pcA <- pd2pc(pdA, pGuess)
  ## Get critical value in one-tailed binomial test:
  crit <- if(is.null(crit))
    findcr(sample.size=ss, alpha=alpha, p0=pGuess, pd0=pd0,
           test=test)
  ## compute power of the test from critical value:
  if(test == "difference") {
    crit <- delimit(crit, lower = 1, upper = ss + 1)
    power <- 1 - pbinom(q=crit - 1, size=ss, prob=pcA)
  }
  else if(test == "similarity") {
    crit <- delimit(crit, lower = 0, upper = ss)
    power <- pbinom(q = crit, size = ss, prob = pcA)
  }
  else ## should never happen
    stop("'test' not recognized")
  return(as.vector(power))
}

discrimPwr <-
  function(pdA, pd0 = 0, sample.size, alpha = 0.05, pGuess = 1/2,
           test = c("difference", "similarity"),
           statistic = c("exact", "normal", "cont.normal"))
### Exported, visible function with argument tests.
{
  ## match and test arguments:
  test <- match.arg(test)
  stat <- match.arg(statistic)
  ss <- sample.size
  stopifnot(is.numeric(pdA) && length(pdA) == 1 &&
            pdA >= 0 && pdA <= 1)
  stopifnot(is.numeric(pd0) && length(pd0) == 1 &&
            pd0 >= 0 && pd0 <= 1)
  stopifnot(is.numeric(ss) && length(ss) == 1 &&
            isTRUE(all.equal(round(ss), ss)) &&
            ss > 0)
  ss <- as.integer(round(ss))
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(is.numeric(pGuess) && length(pGuess) == 1 &&
            pGuess >= 0 && pGuess < 1)
  ## Test admissibility of pdA and pd0:
  if(test == "difference" && pdA < pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA > pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  ## Get appropriate power:
  if(stat == "normal")
    pwr <- normalPwr(pdA = pdA, pd0 = pd0, sample.size = ss,
                     alpha = alpha, pGuess = pGuess, test = test)
  else if(stat == "cont.normal")
    pwr <- normalPwr(pdA=pdA, pd0=pd0, sample.size=ss,
                     alpha=alpha, pGuess=pGuess, test=test,
                     continuity=TRUE)
  else if(stat == "exact")
    pwr <- pdPwr(pdA=pdA, pd0=pd0, sample.size=ss, alpha=alpha,
                 pGuess=pGuess, test=test)
  else stop('"statistic" not recognized')
  return(as.vector(pwr))
}

d.primePwr <-
  function(d.primeA, d.prime0 = 0, sample.size, alpha = 0.05,
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle"),
           test = c("difference", "similarity"),
           statistic = c("exact", "normal", "cont.normal"))
{
  ## Convenience function that simply modifies some arguments and
  ## calls discrimPwr
  newCall <- call <- match.call()
  method <- match.arg(method)
  stopifnot(length(d.primeA) == 1 && is.numeric(d.primeA) &&
            d.primeA >= 0)
  stopifnot(length(d.prime0) == 1 && is.numeric(d.prime0) &&
            d.prime0 >= 0)
  pdA <- coef(rescale(d.prime = d.primeA, method = method))$pd
  pd0 <- coef(rescale(d.prime = d.prime0, method = method))$pd
  newCall$method <- newCall$d.primeA <- newCall$d.prime0 <- NULL
  newCall$pGuess <-
    ifelse(method %in% c("duotrio", "twoAFC"), 1/2, 1/3)
  newCall$pdA <- pdA
  newCall$pd0 <- pd0
  newCall[[1]] <- as.name("discrimPwr")
  return(eval.parent(newCall))
}

