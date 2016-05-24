estimate.model.lnre.zm <- function (model, spc, param.names,
                                     method, cost.function, m.max=15,
                                     debug=FALSE, ...)
{
  if (! inherits(model, "lnre.zm")) stop("argument must be object of class 'lnre.zm'")

  ## root function used to solve implicit equation E[V(N)] = V for B
  B.root.function <- function (B, N, V, model) {
    B <- B
    model <- model$util$update(model, list(B=B))
    E.V <- EV(model, N)
    E.V - V
  }
  
  ## B is computed directly (so that E[V] = V) and not controlled by the optimization procedure
  compute.B <- function(model, spc) {
    alpha <- model$param$alpha
    N <- N(spc)
    V <- V(spc)
    
    ## estimate for the non-exact model can be computed directly
    C <- alpha * V / (Cgamma(1-alpha) * N^alpha)
    B <- ( (1 - alpha) / C ) ^ (1 / (1 - alpha))
    if (B > 1e6) {
      if (debug) warning("Estimate B = ",B," > 1e6 from inexact model, adjusting to 1e6")
      B <- 1e6
    }
    
    ## for the exact model, this estimate leads to E[V] < V (because of the nature of the
    ## approximation), so we have to adjust B appropriately in this case
    if (model$exact) {
      .eps <- 10 * .Machine$double.eps
      upper.B <- B
      lower.B <- upper.B
      .diff <- B.root.function(lower.B, N, V, model)
      if (.diff >= 0) { # if E[V] < V isn't met, don't know how to search for better solution
        if (debug) warning("Can't adjust estimate B = ",B," from inexact model ",
                           "since E[V] - V = ",.diff," > 0")
      }
      else { # now find a lower boundary for the search interval with E[V] > V
        while (lower.B > 0 && .diff < 0) { 
          lower.B <- lower.B / 2
          .diff <- B.root.function(lower.B, N, V, model)
        }
        if (.diff < 0) { # can't find suitable lower boundary -> fall back to inexact estimate
          if (debug) warning("Can't find exact estimate for B, since all B > 0 have ",
                             "E[V] - V <= ",.diff," < 0")
        }
        else {       # now use uniroot() to find B in this range with E[V] = V
          res <- uniroot(B.root.function, lower=lower.B, upper=upper.B, tol=.eps * B,
                         N=N, V=V, model=model)
          if (abs(res$f.root) > .5) {   # estimate for B is not satisfactory
            if (debug) warning("No exact estimate for B found.  Best guess was ",
                               "B = ",res$root," with E[V] - V = ",res$f.root)
          }
          else {     # found satisfactory exact estimate for B -> update value
            B <- res$root
          }
        }
      }
    }
    
    B
  }
  
  if ("B" %in% param.names) {           
    param.names <- param.names[param.names != "B"]
  }
  else {
    warning("parameter B cannot be fixed in ZM estimation, ignoring specified value")
  }

  param.values <- rep(0, length(param.names)) 
  if ("alpha" %in% param.names) {
    tmp <- list(alpha=Vm(spc, 1) / V(spc)) # start with direct ZM estimate, cf. Evert (2004), p. 130
    tmp <- model$util$transform(tmp)
    param.values[param.names == "alpha"] <- tmp$alpha
  }
     
  compute.cost <- function (P.vector, param.names, model, spc, m.max=15, debug=FALSE)
    {
      P.trans <- as.list(P.vector)                     # translate parameter vector into list
      names(P.trans) <- param.names
      P <- model$util$transform(P.trans, inverse=TRUE) # convert parameters to normal scale
      P$B <- 1 # make sure that B is valid (just a dummy, correct value is computed below)
      model <- model$util$update(model, P)             # update model parameters (checks ranges)

      B <- compute.B(model, spc) # compute B from other model parameters and spectrum
      model <- model$util$update(model, list(B=B))
      
      cost <- cost.function(model, spc, m.max)

      if (debug) {
        report <- as.data.frame(model$param)
        report$cost <- round(cost, digits=2)
        rownames(report) <- ""
        print(report)
      }
      cost
    }

  if (length(param.names) > 0) { # at most 1 parameter is estimated -> one-dimensional minimization with NLM
    result <- nlm(compute.cost, param.values, print.level=debug, stepmax=10, steptol=1e-12,
                  param.names=param.names, model=model, spc=spc, m.max=m.max, debug=debug)

    res.code <- result$code
    if (res.code > 3) stop("parameter estimation failed (code ", res.code,")")
    if (res.code == 3) warning("estimated parameter values may be incorrect (code 3)")

    P.estimate <- as.list(result$estimate)
    names(P.estimate) <- param.names

    model <- model$util$update(model, P.estimate, transformed=TRUE)
  }
  ## if B was the only parameter to be estimated,  we can just compute it directly

  B <- compute.B(model, spc)
  model <- model$util$update(model, list(B=B))
  
  model$gof <- lnre.goodness.of.fit(model, spc, n.estimated=length(param.names) + 1)
  
  model
}

