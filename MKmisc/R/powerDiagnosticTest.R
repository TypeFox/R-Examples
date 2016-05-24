power.diagnostic.test <- function(sens = NULL, spec = NULL,
                                  n = NULL, delta = NULL, sig.level = 0.05,
                                  power = NULL, prev = NULL, 
                                  method = c("exact", "asymptotic"),
                                  NMAX = 1e4){
  if(sum(sapply(list(sens, spec), is.null)) != 1) 
    stop("exactly one of 'sens', and 'spec' must be NULL")
  if(sum(sapply(list(n, delta, sig.level, power), is.null)) != 1) 
    stop("exactly one of 'n', 'delta', 'sig.level', and 'power' must be NULL")
  if(!is.null(delta) && !is.numeric(delta) || any(0 > delta | delta > 1))
    stop("'delta' must be numeric in [0, 1]")
  if(!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  if(!is.null(power) && !is.numeric(power) || any(0 > power | power > 1))
    stop("'power' must be numeric in [0, 1]")
  if(!is.null(prev) && !is.numeric(prev) || any(0 > prev | prev > 1))
    stop("'prev' must be numeric in [0, 1]")
  
  if(is.null(spec))
    prob <- sens
  else
    prob <- spec
  if(is.null(prev)) 
    pr <- 0.5
  else
    pr <- prev
  
  if(!is.null(delta)){
    prob.delta <- prob - delta
    if(!is.null(prob.delta) && !is.numeric(prob.delta) || any(0 > prob.delta | prob.delta > 1))
      stop("'sens'-'delta' resp. 'spec'-'delta' must be numeric in [0, 1]")    
  }
  method <- match.arg(method)
  method.nr <- switch(method, exact = 1, asymptotic = 2)
  
  ## exact
  if(method.nr == 1){
    if(is.null(n)){
      ns <- 2:NMAX
      lows <- qbinom(1-power, size = ns, prob = prob) - 1
      crit <- pbinom(lows, size = ns, prob = prob.delta, lower.tail = FALSE) - sig.level
      ind <- max(which(crit > 0)) + 1
      n <- ns[ind]
    }else if(is.null(delta)){
      fun.delta <- function(delta){
        lo <- qbinom(1-power, size = n, prob = prob) - 1
        pbinom(lo, size = n, prob = prob-delta, lower.tail = FALSE) - sig.level
      }
      delta <- uniroot(f = fun.delta, interval = c(1e-10,prob-1e-10))$root
    }else if(is.null(power)){
      fun.power <- function(power){
        lo <- qbinom(1-power, size = n, prob = prob) - 1
        pbinom(lo, size = n, prob = prob.delta, lower.tail = FALSE) - sig.level
      }
      power <- uniroot(f = fun.power, interval = c(1e-10,1-1e-10))$root      
    }else if(is.null(sig.level)){
      fun.sig.level <- function(sig.level){
        lo <- qbinom(1-power, size = n, prob = prob) - 1
        pbinom(lo, size = n, prob = prob.delta, lower.tail = FALSE) - sig.level
      }
      sig.level <- uniroot(f = fun.sig.level, interval = c(1e-10,1-1e-10))$root            
    }else stop("internal error")
  }
  
  ## asymptotic
  if(method.nr == 2){
    if(is.null(n))
      n <- (qnorm(power)*sqrt(prob*(1-prob)) 
            + qnorm(1-sig.level)*sqrt(prob.delta*(1-prob.delta)))^2/delta^2
    else if(is.null(delta)){
      fun.delta <- function(delta){ 
        sqrt(n)*delta - (qnorm(power)*sqrt(prob*(1-prob)) 
                         + qnorm(1-sig.level)*sqrt((prob - delta)*(1 - prob + delta)))
      }
      delta <- uniroot(f = fun.delta, interval = c(1e-10,prob-1e-10))$root
    }else if(is.null(power))
      power <- pnorm((sqrt(n)*delta - qnorm(1-sig.level)*sqrt(prob.delta*(1-prob.delta)))/sqrt(prob*(1-prob)))
    else if(is.null(sig.level))
      sig.level <- 1 - pnorm((sqrt(n)*delta - qnorm(power)*sqrt(prob*(1-prob)))/sqrt(prob.delta*(1-prob.delta)))
    else stop("internal error")
  }
  
  if(is.null(spec))
    n1 <- n*(1-pr)/pr
  else
    n1 <- n*pr/(1-pr)        
  METHOD <- paste("Diagnostic test", switch(method, exact = "exact", asymptotic = "asymptotic"), 
                  "power calculation")
  if(is.null(spec)){
    NOTE <- "n is number of cases, n1 is number of controls"
    res <- structure(list(sens = sens, n = n, n1 = n1, delta = delta, 
                          sig.level = sig.level, power = power, prev = prev, 
                          note = NOTE, method = METHOD), 
            class = "power.htest")  
  }else{
    NOTE <- "n is number of controls, n1 is number of cases"
    res <- structure(list(sens = sens, n = n, n1 = n1, delta = delta, 
                          sig.level = sig.level, power = power, prev = prev, 
                          note = NOTE, method = METHOD), 
                     class = "power.htest")
  }
  res
}
