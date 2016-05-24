normalSS <-
  function(pdA, pd0 = 0, target.power = 0.9, alpha = 0.05,
           pGuess = 1/2, test = c("difference", "similarity"),
           type=c("normal", "cont.normal", "lwr", "upr", "upr.EJ",
             "user"), c = 1/2)
{
### Non-user-visible function assuming admissible arguments.
### EJ.type: compute the bound from Ennis and Jesionka (2011) or use
### my own version.

  adjust.c <- function(alpha) {
    ## Adjust the c-constant as a function of alpha (or beta) when
    ## computing lwr and upr.
    if(alpha >= 10^-2) return(-1)
    -1 - 1 * ceiling(-log10(alpha * 100))
  }
  test <- match.arg(test)
  type <- match.arg(type)
  ## Cannot compute normal ss if (alpha, beta) >= 0.5 when c!=0.
  ## We allow this using type="user", though.
  ## Also note that type="normal" is ok.
  if((target.power <= 0.5 || alpha >= 0.5) &&
     type %in% c("cont.normal", "lwr", "upr", "upr.EJ"))
    stop("Normal approximation requires target.power > 0.5 and alpha < 0.5")  
  pow <- target.power
  pc0 <- pd2pc(pd=pd0, Pguess=pGuess)
  pcA <- pd2pc(pd=pdA, Pguess=pGuess)
  if(test == "similarity") {
    ## Frame similarity test as a difference test:
    ## Swap (pcA, pc0) and (alpha, beta):
    pcA.temp <- pcA
    pcA <- pc0
    pc0 <- pcA.temp
    alpha.temp <- alpha
    alpha <- 1-pow
    pow <- 1-alpha.temp
  }
  if(type == "upr.EJ") alpha <- alpha/3
  ## Solve 'pow(n) = target' for n (2nd degree polynomial):
    a <- pc0 - pcA
    b <- qnorm(1 - alpha) * sqrt(pc0 * (1 - pc0)) -
      qnorm(1 - pow) * sqrt(pcA * (1 - pcA))
    c <- switch(type, "normal" = 0,
                "cont.normal" = .5,
                "lwr" = adjust.c(alpha),
                "upr" = 1 - adjust.c(1-pow),
                "upr.EJ" = .5,
                "user" = c)
    D <- max(0, b^2 - 4 * a * c)
    x <- (-b - sqrt(D)) / (2 * a)
    n <- max(1, round(x^2))
  return(as.vector(n))
}    

pdSS <-
  function(pdA, pd0 = 0, target.power = 0.90, alpha = 0.05,
           pGuess = 1/2, test = c("difference", "similarity"),
           statistic = c("exact", "stable.exact", "both.exact"))
{
  test <- match.arg(test)
  stat <- match.arg(statistic)
  pow <- target.power
  pcA <- pd2pc(pdA, pGuess)
  sample.size <- double(length=0L)
  ## Upper and lower bounds based on the normal approximations: 
  n.lwr <- normalSS(pdA = pdA, pd0 = pd0, target.power = target.power,
                    alpha = alpha, pGuess = pGuess, test = test,
                    type="lwr") ## -5
  n.upr <- normalSS(pdA = pdA, pd0 = pd0, target.power = target.power,
                    alpha = alpha, pGuess = pGuess, test = test,
                    type="upr") ## +5
  ## Power at these bounds:
  lwr.pow <- pdPwr(pcA=pcA, pd0=pd0, sample.size=n.lwr, alpha=alpha,
                   pGuess=pGuess, test=test)
  upr.pow <- pdPwr(pcA=pcA, pd0=pd0, sample.size=n.upr, alpha=alpha,
                   pGuess=pGuess, test=test)
  ## Test that target.power is between lwr.pow and upr.pow:
  if(lwr.pow >= pow) stop("Unable to determine lower bound on n")
  if(upr.pow <= pow) stop("Unable to determine upper bound on n")
  ## (First) exact and stable exact sample sizes:
  if(stat %in% c("exact", "both.exact")) {
    ## Compute "first" exact sample size:
    ## Step up from n.lwr until power > target.power:
    for(i in (n.lwr + 1):n.upr) {
      lwr.pow <- pdPwr(pcA=pcA, pd0=pd0, sample.size=i, alpha=alpha,
                       pGuess=pGuess, test=test)
      if(lwr.pow > pow) break
    }
    sample.size <- c(sample.size, structure(i, names="exact"))
  }
  if(stat %in% c("stable.exact", "both.exact")) {
    ## Compute stable exact sample size:
    ## Step down from n.upr until power < target.power:
    for(i in (n.upr - 1):n.lwr) {
      upr.pow <- pdPwr(pcA=pcA, pd0=pd0, sample.size=i, alpha=alpha,
                       pGuess=pGuess, test=test)
      if(upr.pow < pow) break
    }
    sample.size <- c(sample.size,
                     structure(i + 1, names="stable.exact"))
  }
  return(as.vector(sample.size))
}

discrimSS <-
  function(pdA, pd0 = 0, target.power = 0.90, alpha = 0.05,
           pGuess = 1/2, test = c("difference", "similarity"),
           statistic = c("exact", "stable.exact", "both.exact",
             "normal", "cont.normal")) 
{
  ## Match and test arguments:
  test <- match.arg(test)
  call <- match.call()
  stat <- match.arg(statistic)
  stopifnot(is.numeric(pdA) && length(pdA) == 1 &&
            pdA >= 0 && pdA <= 1)
  stopifnot(is.numeric(pd0) && length(pd0) == 1 &&
            pd0 >= 0 && pd0 <= 1)
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(is.numeric(target.power) && length(target.power) == 1 && 
            target.power > 0 && target.power < 1)
  stopifnot(is.numeric(pGuess) && length(pGuess) == 1 &&
            pGuess >= 0 && pGuess < 1)
  if(test == "difference" && pdA <= pd0)
    stop("pdA has to be larger than pd0 for difference tests")
  if(test == "similarity" && pdA >= pd0)
    stop("pdA has to be less than pd0 for similarity tests")
  ## Normal approximation to sample size:
  Stat <- ifelse(stat == "normal", stat, "cont.normal")
  n <- normalSS(pdA = pdA, pd0 = pd0, target.power = target.power,
                alpha = alpha, pGuess = pGuess, test = test,
                type=Stat)
  ## Return appropriate sample size:
  if(stat %in% c("normal", "cont.normal")) return(n)
  if(n > 1e5) {
    warning(paste("sample size probably > 1e4 and 'exact' option is",
                  "not available\n",
                  "using normal approximation instead"))
    return(as.vector(n))
  }
  n <- pdSS(pdA = pdA, pd0 = pd0, target.power = target.power,
            alpha = alpha, pGuess = pGuess, test = test,
            statistic=stat)
  return(as.vector(n))
}

d.primeSS <- 
  function(d.primeA, d.prime0 = 0, target.power = 0.90, alpha = 0.05,
           method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle"), 
           test = c("difference", "similarity"),
           statistic = c("exact", "stable.exact", "both.exact",
             "normal", "cont.normal"))
{
  ## Convenience function that simply modifies some arguments and
  ## calls discrimSS
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
  newCall[[1]] <- as.name("discrimSS")
  return(eval.parent(newCall))
}
