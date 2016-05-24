# Implements MM algorithm with generalized gradient descent (and Nesterov's
# method) using ADMM to solve prox function.

spcov <- function(Sigma, S, lambda, step.size, nesterov=TRUE,
               n.outer.steps=1e4, n.inner.steps=1e4,
               tol.outer=1e-4, thr.inner=1e-2,
               backtracking=0.2, trace=0) {
  # Performs MM to optimize the non-convex objective function
  # Args:
  #  Sigma: an initial guess for Sigma
  #  lambda: either a scalar or a matrix of the same dimensions as Sigma
  #  nesterov: indicates whether to use Nesterov or standard generalized gradient
  #                    descent to perform inner loop optimization
  #  tol.outer: convergence threshold for outer (MM) loop.
  #  thr.inner: convergence threshold for inner (GG/Nesterov) loop.
  #             stops when mean absolute change in Sigma is
  #             less than thr.inner * mean(abs(S))
  #  backtracking: see "GGDescent"
  if (all(lambda == 0)) {
    cat("Skipping MM.  Solution is S!", fill=T)
    return(list(n.iter=0, Sigma=S, obj=ComputeObjective(S, S, lambda)))
  }
  stopifnot(lambda >= 0)
  if (trace > 0) {
    cat("---", fill=T)
    cat(ifelse(nesterov, "using Nesterov, ", ""))
    cat(ifelse(backtracking, "backtracking line search", ""), fill=T)
    cat("---", fill=T)
  }
  mean.abs.S <- mean(abs(S))
  if (min(eigen(Sigma, symmetric=T, only.values=T)$val) < 1e-5)
    warning("Starting value is nearly singular.")

  del <- ComputeDelta(S, lambda, trace=trace-1) # get a lower bound on minimum eval

  objective <- ComputeObjective(Sigma, S, lambda)
  if (trace > 0)
    cat("objective: ", objective, fill=T)
  n.iter <- NULL # number of inner iterations on each step
  for (i in seq(n.outer.steps)) {
    Sigma0 <- Sigma
    if (trace > 0)
      cat("step size given to GGDescent/Nesterov:", step.size, fill=T)
    gg <- GGDescent(Sigma=Sigma, Sigma0=Sigma0, S=S, lambda=lambda,
                    del=del, nsteps=n.inner.steps,
                    step.size=step.size,
                    nesterov=nesterov,
                    tol=thr.inner * mean.abs.S,
                    trace=trace - 1,
                    backtracking=backtracking)
    Sigma <- gg$Sigma
    objective <- c(objective, ComputeObjective(Sigma, S, lambda))
    if (trace > 0) {
      cat("objective: ", objective[length(objective)],
          " (", gg$niter, "iterations, max step size:",
          max(gg$step.sizes), ")",
          fill=T)
    }
    if (backtracking) {
      if (max(gg$step.sizes) < step.size * backtracking ^ 2) {
        step.size <- step.size * backtracking
        if (trace > 0)
          cat("Reducing step size to", step.size, fill=T)
      }
    }
    n.iter <- c(n.iter, gg$niter)
    if(objective[i + 1] > objective[i] - tol.outer) {
       cat("MM converged in", i, "steps!", fill=T)
       break
    }
  }
  list(n.iter=n.iter, Sigma=gg$Sigma, obj=objective)
}

GGDescent <- function(Sigma, Sigma0, S, lambda, del, nsteps,
                                     step.size, nesterov=FALSE, backtracking=FALSE,
                                     tol = 1e-3, trace=0) {
  # solves the problem
  # Min_{Sigma} tr(solve(Sigma0,Sigma)) + tr(solve(Sigma,S))
  #                                     + ||lambda * Sigma||_1
  # using Nesterov's method with backtracking.
  # Note: We wish to solve this with the constraint that Sigma pd.
  # However, this algorithm does not impose this constraint.
  # Args:
  #  Sigma: an initial guess for Sigma
  #  Sigma0, S, lambda, del: parameters of the optimization problem
  #  nsteps: number of generalized gradient steps to take
  #  nesterov: TRUE/FALSE, indicates whether to take Nesterov vs. standard gen grad steps.
  #  backtracking: if FALSE, then fixed step size used.  If numeric and in
  #                (0,1), this is the beta parameter of backtracking.
  #                Usually, beta is in (0.1, 0.8).
  #  tol:  convergence threshold.  Stops when mean(abs(Sigma-Sigma.last)) < tol
  if (backtracking) {
    beta <- backtracking
    if (beta <= 0 | beta >= 1)
      stop("Backtracking parameter beta must be in (0,1).")
  }
  tt <- step.size
  converged <- FALSE
  exit <- FALSE
  obj.starting <- ComputeObjective(Sigma, S, lambda)
  Sigma.starting <- Sigma
  Omega <- Sigma
  Sigma.last <- Sigma
  ttts <- NULL
  ttt <- tt # note: as in Beck & Teboulle, step size only gets smaller.
  for (i in seq(nsteps)) {
    inv.Sigma0 <- solve(Sigma0)
    log.det.Sigma0 <- LogDet(Sigma0)
    grad.g <- ComputeGradientOfg(Omega, S, Sigma0, inv.Sigma0=inv.Sigma0)
    grad.g <- (grad.g + t(grad.g)) / 2 # make sure this stays symmetric
    g.omega <- g(Omega, S, Sigma0,
                 inv.Sigma0=inv.Sigma0, log.det.Sigma0=log.det.Sigma0)
    # backtracking line search:
    while (backtracking) {
      #soft.thresh <- SoftThreshold(Omega - ttt * grad.g, lambda * ttt)
      soft.thresh <- ProxADMM(Omega - ttt * grad.g, del, 1, P=lambda*ttt, rho=.1)$X
      gen.grad.g <- (Omega - soft.thresh) / ttt
      left <- g(soft.thresh, S, Sigma0,
                inv.Sigma0=inv.Sigma0, log.det.Sigma0=log.det.Sigma0)
      right <- g.omega - ttt * sum(grad.g * gen.grad.g) + ttt * sum(gen.grad.g ^ 2) / 2
      if (is.na(left) || is.na(right)) {
        print("left or right is NA.")
        browser()
      }
      if (left <= right) {
        # accept this step size
        Sigma <- soft.thresh
        ttts <- c(ttts, ttt)
        # check for convergence
        if (mean(abs(Sigma - Sigma.last)) < tol) {
          converged <- TRUE
          break # note: this break only ends the backtracking loop
        }
        if (nesterov)
          Omega <- Sigma + (i - 1) / (i + 2) * (Sigma - Sigma.last)
        else
          Omega <- Sigma
        Sigma.last <- Sigma
        if (trace > 0)
          cat("--true objective:", ComputeObjective(Sigma, S, lambda), fill=T)
        if (trace > 0)
          cat(i, ttt, " ")
        break
      }
      ttt <- beta * ttt
      if (ttt < 1e-15) {
        cat("Step size too small: no step taken", fill=T)
        exit <- TRUE
        break
      }
    }
    if (!backtracking) {
      #Sigma <- SoftThreshold(Sigma - ttt * grad.g, lambda * ttt)
      Sigma <- ProxADMM(Sigma - ttt * grad.g, del, 1, P=lambda*ttt, rho=.1)$X
      # check for convergence:
      if (mean(abs(Sigma - Sigma.last)) < tol)
        converged <- TRUE
      if (nesterov)
        Omega <- Sigma + (i - 1)/(i + 2) * (Sigma - Sigma.last)
      else
        Omega <- Sigma
      Sigma.last <- Sigma
    }
    if (converged) {
      if(trace > 0) {
        cat("--GG converged in", i, "steps!")
        if (backtracking)
          cat(" (last step size:", ttt, ")", fill=T)
        else
          cat(fill=T)
      }
      break
    }
    if (exit) {
      break
    }
  }
  obj.end <- ComputeObjective(Sigma, S, lambda)
  if (obj.starting < obj.end) {
    if (nesterov) {
      cat("Objective rose with Nesterov.  Using generalized gradient instead.", fill=T)
      return(GGDescent(Sigma=Sigma.starting, Sigma0=Sigma0, S=S, lambda=lambda,
                       del=del, nsteps=nsteps, step.size=step.size,
                       nesterov=FALSE,
                       backtracking=backtracking,
                       tol = tol, trace=trace))
    }
    
    browser()
    cat("--Returning initial Sigma since GGDescent/Nesterov did not decrease objective", fill=T)
    Sigma <- Sigma.starting
  }

  list(Sigma=Sigma, niter=i, step.sizes=ttts)
}


ComputeObjective <- function(Sigma, S, lambda) {
  # the original non-convex problem's objective function
  -2 * ComputeLikelihood(Sigma, S) + ComputePenalty(Sigma, lambda)
}
