### INTERNAL DATA OBJECTS
# data    list with the following components:
#           X    n-by-p matrix
#           y    n-by-q matrix
#           XtX  t(X) %*% X
#           Xty  t(X) %*% y
# state   current state of the MCMC samples; list with the following components:
#           indic.grp    G-by-1 logical matrix if posterior uses indic.grp;
#                        otherwise NULL or missing
#           indic.var        p-by-q logical matrix
#           tau          scalar
#           S2           S_Z^2, etc
#           OmegaInvDet  log of determinant of Omega_Z^{-1} etc
#                        if NOT using g-prior; otherwise NULL or missing
#           XDet         log of determinant of X_Z, etc, if using
#                        determinantal prior; otherwise NULL or missing
# TRANSLATOR between code and the notation in Stell and Sabatti (2015):
#   CODE      PAPER
#   indic.var Z
#   indic.grp W for Across Traits, G for Across Sites; throughout the code, grp
#             stands for combining information across either traits or sites
#   omega     omega for basic, nu for other priors
#   omega.grp omega_W or omega_G

ptycho <- function(X, y, initStates, groups=NULL,
                  tau.min=0.01, tau.max=10, tau.sd=(tau.max-tau.min)/4,
                  doGPrior=TRUE, doDetPrior=FALSE, prob.varadd=0.5,
                  isOmegaFixed=FALSE, omega=NULL, omega.grp=NULL,
                  probs.grp=NULL, rho.alpha=10, rho.lambda=rho.alpha,
                  only.means=FALSE, nburn=0, nthin=1, nSavePerChain,
                  parallel.chains=FALSE, random.seed=NULL) {
  doGrpIndicator <- !is.null(initStates[[1]]$indic.grp)
  if (isOmegaFixed) {
    if (is.null(omega)) {
      stop("ptycho: omega must be specified when isOmegaFixed is TRUE")
    } else if (doGrpIndicator && is.null(omega.grp)) {
      stop("ptycho: Using indic.grp and isOmegaFixed TRUE requires omega.grp")
    }
  }
  data <- list(X=X, y=y, XtX=t(X)%*%X, Xty=t(X)%*%y)
  if (!doGrpIndicator) {
    if (!is.null(omega.grp)) {
      stop("ptycho: omega.grp given but no indic.grp in initial states")
    }
    if (is.null(omega)) {
      x <- if (ncol(y)==1) 1 else ncol(X)
      omega <- cbind(A=rep(1,x),B=rep(1,x))
    }
  } else {
    G <- if (is.null(groups)) ncol(X) else length(groups$group2var)
    if (is.null(probs.grp)) probs.grp <- c(0.25,0.5,0.25)
    if (is.null(omega)) omega <- cbind(A=rep(1,G),B=rep(1,G))
    if (is.null(omega.grp)) omega.grp <- c(A=1,B=1)
  }
  if (is.logical(only.means)) {
    only.means <- if (only.means) seq_len(nSavePerChain) else NULL
  } else {
    nSavePerChain <- max(only.means)
  }
  params <- list(nSavePerChain=nSavePerChain, nburn=nburn, nthin=nthin,
                 only.means=only.means,
                 prob.varadd=prob.varadd, probs.grp=probs.grp, groups=groups,
                 isOmegaFixed=isOmegaFixed, omega=omega, omega.grp=omega.grp,
                 rho.alpha=rho.alpha, rho.lambda=rho.lambda,
                 tau.min=tau.min, tau.max=tau.max, tau.sd=tau.sd)
  st <- llply(initStates,
              function(x) { initAuxVars(data, x, doGPrior, doDetPrior) })
  chainIterator <- (if (parallel.chains)
                      if (is.null(random.seed)) "chainLoopParallel"
                      else "chainLoopRNG"
                    else "chainLoop")
  z <- do.call(chainIterator,
               list(data=data, params=params, initStates=st,
                    random.seed=random.seed))
  if (is.null(only.means)) {
    z <- mcmc.list(z)
  } else {
    z <- cbind(chain=rep(seq_along(z), each=length(params$only.means)), 
               do.call("rbind", z))
    class(z) <- "ptycho"
  }
  attr(z,"params") <- params
  z
}

chainLoop <- function(data, params, initStates, ncpu=NULL, random.seed=NULL) {
  nchains <- length(initStates)
  z <- vector("list", nchains)
  if (!is.null(random.seed)) set.seed(random.seed)
  for (nn in seq_len(nchains)) {
    state <- initStates[[nn]]
    z[[nn]] <- mcmcLoop(data, params, state)
  }
  z
}

# The following is in batch.R and applies here as well.
#utils::globalVariables(c("nn"))

# NOT REPRODUCIBLE.
# random.seed is ignored but included to make do.call() in ptycho() simple.
chainLoopParallel <- function(data, params, initStates, ncpu, random.seed=NULL) {
  nchains <- length(initStates)
  checkParallel()
  z <- foreach::"%dopar%"(foreach::foreach(nn=seq_len(nchains)),
                         {
                           state <- initStates[[nn]]
                           mcmcLoop(data, params, state)
                         })
  z
}

# REPRODUCIBLE requires doRNG
chainLoopRNG <- function(data, params, initStates, ncpu, random.seed) {
  nchains <- length(initStates)
  checkParallel(checkRNG=TRUE)
  set.seed(random.seed)
  z <- doRNG::"%dorng%"(foreach::foreach(nn=seq_len(nchains)),
                         {
                           state <- initStates[[nn]]
                           mcmcLoop(data, params, state)
                         })
  z
}

# Finish filling in the initial states
initAuxVars <- function(data, state, doGPrior, doDet) {
  q <- ncol(state$indic.var)
  state$S2 <- rep(NA,q)
  if (!doGPrior) state$OmegaInvDet <- rep(NA,q)
  if (doDet) state$XDet <- rep(NA,q)
  for (k in seq_len(q)) state <- updateAuxVars(data, state, k)
  state
}

# Update S2[k], OmegaInvDet[k], and XDet[k].
updateAuxVars <- function(data, state, k) {
  state$S2[k] <- computeSg2(data, state, k)
  if (!is.null(state$OmegaInvDet)) {
    state$OmegaInvDet[k] <- computeLogDet(data, state, k, useTau=TRUE)
  }
  if (!is.null(state$XDet)) {
    state$XDet[k] <- computeLogDet(data, state, k, useTau=FALSE)
  }
  state
}

computeSg2 <- function(data, state, k) {
  yty <- sum(data$y[,k]^2)
  indic.var <- state$indic.var[,k]
  if (!any(indic.var)) return(yty)
  OmegaInv <- data$XtX[indic.var,indic.var,drop=FALSE]
  if (is.null(state$OmegaInvDet)) {
    n <- nrow(data$X)
    OmegaInv <- (1+n*state$tau^2) * OmegaInv / (n*state$tau^2)
  } else {
    OmegaInv <- OmegaInv + diag(nrow=nrow(OmegaInv))/state$tau^2
  }
  S2 <- data$Xty[indic.var,k] %*% solve(OmegaInv, data$Xty[indic.var,k])
  S2 <- yty - S2
}

computeLogDet <- function(data, state, k, useTau) {
  indic.var <- state$indic.var[,k]
  svar <- sum(indic.var)
  if (svar == 0) return(0)
  z <- data$XtX[indic.var,indic.var]
  z <- if (useTau) z + diag(nrow=svar)/state$tau^2 else z/nrow(data$XtX)
  z <- eigen(z, symmetric=TRUE, only.values=TRUE)
  z <- sum(log(z$values))
}

### MCMC LOOP
mcmcLoop <- function(data, params, state) {
  z <- initMCMCArray(state, params$nSavePerChain, params$only.means,
                     colnames(data$y))
  if (!is.null(params$only.means)) zsum <- rep(0, ncol(z)-1)
  nrow <- 1
  for (nsmpl in seq_len(params$nthin*params$nSavePerChain+params$nburn)) {
    state <- mhIndicators(data, params, state)
    state <- mhTau(data, params, state)
    if (nsmpl > params$nburn) {
      nsv <- nsmpl - params$nburn
      if (is.null(params$only.means)) {
        if (nsv %% params$nthin == 0) {
          zrow <- c(as.numeric(state$indic.grp), as.numeric(state$indic.var),
                    state$tau)
          z[nrow,] <- zrow
          nrow <- nrow + 1
        }
      } else {
        t <- state$tau
        zrow <- c(as.numeric(state$indic.grp), as.numeric(state$indic.var),
                  t, t^2)
        zsum <- zsum + zrow
        if (nsv %in% params$only.means) {
          z[nrow,] <- c(nsv, zsum/nsv)
          nrow <- nrow + 1
        }
      }
    }
  }
  if (is.null(params$only.means)) z <- mcmc(z, start=1, thin=1)
  z
}

### INITIALIZE ARRAY FILLED IN BY mcmcLoop()
initMCMCArray <- function(state, nrow, only.means, ynames=NULL) {
  p <- nrow(state$indic.var)
  q <- ncol(state$indic.var)
  if (is.null(ynames)) ynames <- createYNames(q)
  stnames <- (if (is.null(state$indic.grp)) NULL
              else if (q==1) paste("indic.grp[", seq_along(state$indic.grp),
                                   ",", ynames, "]", sep="")
              else paste("indic.grp[", seq_along(state$indic.grp), "]",
                         sep=""))
  stnames <- c(stnames,
               paste("indic.var[", rep(seq_len(p),times=q), ",",
                     rep(ynames,each=p), "]", sep=""))
  stnames <- c(stnames, "tau")
  if (is.null(only.means)) {
    z <- matrix(nrow=nrow, ncol=length(stnames))
    colnames(z) <- stnames
  } else {
    z <- matrix(nrow=length(only.means), ncol=2+length(stnames))
    colnames(z) <- c("iter", stnames, "tau2")
  }
  z
}

### APPLY METROPOLIS-HASTINGS TO indic.var and indic.grp
mhIndicators <- function(data, params, state) {
  # Only loop over columns when using multiple traits but not indic.grp
  nloop <- if (is.null(state$indic.grp)) ncol(state$indic.var) else 1
  for (k in seq_len(nloop)) {
    # Compute proposals for indic.grp and indic.var
    if (is.null(state$indic.grp)) {
      grpprop <- NULL
      varprop <- drawIndicVar(state$indic.var[,k], params$prob.varadd,
                              if (params$isOmegaFixed) params$omega[,k]
                              else NULL)
      indic.var <- state$indic.var
      indic.var[,k] <- varprop$proposal
      varprop$proposal <- indic.var
    } else {
      grpprop <- drawIndicGrp(state$indic.grp, params$probs.grp, params$groups)
      if (grpprop$dir == 0) {
        varprop <- drawIndicVarForGrpUnchanged(state$indic.var, grpprop$ind,
                                               params)
      } else if (grpprop$dir == 1) {
        varprop <- drawIndicVarForGrpAdded(state$indic.var, grpprop$ind, params)
      } else {
        varprop <- setIndicVarForGrpDrop(state$indic.var, grpprop$ind,
                                         params$groups)
      }
    }
    # Update state.
    stprop <- updateStateForIndicators(data, state, grpprop, varprop, k)
    # The probability that the proposal is accepted is min(1,r).
    r <- rIndicators(params, nrow(data$X), state, stprop, grpprop, varprop, k)
    if (didChange <- rbinom(1,1,min(1,r))) state <- stprop
  }
  state
}

### APPLY METROPOLIS-HASTINGS TO tau
mhTau <- function(data, params, state) {
  tauprop <- -Inf
  while (tauprop <= params$tau.min || tauprop >= params$tau.max) {
    tauprop <- rnorm(1, mean=state$tau, sd=params$tau.sd)
  }
  # Update state.
  stprop <- updateStateForTau(data, state, tauprop)
  # The probability that the proposal is accepted is min(1,r).
  r <- rTau(params, nrow(data$X), state, stprop)
  if (didChange <- rbinom(1,1,min(1,r))) state <- stprop
  state
}

# Compute contribution to r from S2, OmegaInvDet, and XDet.
# This is used by both rIndicators() and rTau().
computeRAuxVars <- function(params, n, stold, stnew, k) {
  rln <- log(2 * params$rho.lambda + stold$S2[k])
  rln <- rln - log(2 * params$rho.lambda + stnew$S2[k])
  rln <- (0.5*n + params$rho.alpha) * rln
  if (!is.null(stnew$OmegaInvDet)) {
    rln <- rln + 0.5 * (stold$OmegaInvDet[k] - stnew$OmegaInvDet[k])
  }
  # This is pointless in MH for tau, but I don't expect to be using
  # determinantal priors.
  if (!is.null(stnew$XDet)) rln <- rln + stnew$XDet[k] - stold$XDet[k]
  rln
}

### FUNCTIONS FOR INDICATOR VARIABLES
# Input x is only the relevant part of indic.var.  With probability prob, add
# one covariate to model; otherwise remove one.
drawIndicVar <- function(x, prob, omega) {
  xabs <- sum(x); xlen <- length(x)
  if (xabs == 0) {
    jmpdir <- 1
    jmpposs <- seq_len(xlen)
    rfwd <- 1/xlen
    rrev <- if (xlen==1) 1 else 1-prob
  } else if (xabs == xlen) {
    jmpdir <- -1
    jmpposs <- seq_len(xlen)
    rfwd <- 1/xlen
    rrev <- if (xlen==1) 1 else prob
  } else {
    jmpdir <- rbinom(1,1,prob)
    jmpposs <- which(x==(1-jmpdir))
    jmpdir <- 2*jmpdir - 1
    if (jmpdir == 1) {
      rfwd <- prob / (xlen - xabs)
      rrev <- if (xabs == xlen-1) 1/xlen else (1-prob)/(xabs+1)
    } else {
      rfwd <- (1-prob) / xabs
      rrev <- if (xabs == 1) 1/xlen else prob/(xlen-xabs+1)
    }
  }
  # sample() has a "feature" that it chooses from 1:x if length(x)==1.
  if (length(jmpposs) == 1) {
    ind <- jmpposs
  } else {
    ind <- sample(jmpposs, 1)
  }
  x[ind] <- !x[ind]
  list(proposal=x, dir=jmpdir, ind=ind, rfwd=rfwd, rrev=rrev, omega=omega[ind])
}

drawIndicVarForGrpUnchanged <- function(indic.var, grpind, params) {
  if (is.null(params$groups)) {
    vargrp <- indic.var[grpind,]
    omgrp <- params$omega[grpind,]
  } else {
    inds <- params$groups$group2var[[grpind]]
    vargrp <- indic.var[inds]
    omgrp <- params$omega[inds]
  }
  varprop <- drawIndicVar(vargrp, params$prob.varadd, omgrp)
  if (is.null(params$groups)) {
    indic.var[grpind,] <- varprop$proposal
  } else {
    indic.var[inds] <- varprop$proposal
  }
  varprop$proposal <- indic.var
  varprop$omega <- omgrp
  varprop
}

drawIndicVarForGrpAdded <- function(indic.var, grpind, params) {
  if (params$isOmegaFixed) {
    if (is.null(params$groups)) {
      probs <- params$omega[grpind,]
    } else {
      inds <- params$groups$group2var[[grpind]]
      probs <- params$omega[inds]
    }
    indsel <- sapply(probs, function(x) { rbinom(1,1,prob=x) })
  } else {
    if (is.null(params$groups)) {
      grpsize <- ncol(indic.var)
    } else {
      inds <- params$groups$group2var[[grpind]]
      grpsize <- length(inds)
    }
    if (grpsize == 1) {
      indsel <- 1
    } else {
      z <- 0:grpsize
      probs <- beta(params$omega[grpind,"A"] + z,
                    params$omega[grpind,"B"] + grpsize - z)
      probs <- probs * choose(grpsize, z)
      nsel <- sample(z, size=1, prob=probs)
      indsel <- sample.int(grpsize, size=nsel)
    }
  }
  if (is.null(params$groups)) {
    indic.var[grpind,indsel] <- TRUE
  } else {
    indic.var[inds[indsel]] <- TRUE
  }
  list(proposal=indic.var, dir=1, ind=indsel)
}

# No choice here
setIndicVarForGrpDrop <- function(indic.var, grpind, groups) {
  if (is.null(groups)) {
    indsel <- which(indic.var[grpind,])
    indic.var[grpind,indsel] <- FALSE
  } else {
    inds <- groups$group2var[[grpind]]
    indsel <- which(indic.var[inds])
    indic.var[inds[indsel]] <- FALSE
  }
  list(proposal=indic.var, dir=-1, ind=indsel)
}

# Draw indic.grp.  Add or remove a group or leave unchanged.
# When I changed my variable names so that they didn't conflict so egregiously
# with the paper, I decided to leave Z0 and Zh alone inside this function.
drawIndicGrp <- function(indic.grp, probs, groups) {
  sgrp <- sum(indic.grp); ngrp <- length(indic.grp)
  # Zh[i] is TRUE for those entries in indic.grp that we'll choose from if
  # jmpdir==0 (Zh stands for Z-hold)
  if (is.null(groups)) {
    Zh <- indic.grp
    Z0 <- sgrp
  } else {
    Zh <- (indic.grp & (groups$sizes > 1))
    Z0 <- sum(Zh)
  }
  # If jmpdir==0, I only need these values in order to populate returned list.
  rfwd <- NA; rrev <- NA
  if (sgrp == 0) {
    jmpdir <- 1
  } else if (sgrp == ngrp) {
    # Assume some group has more than one member
    jmpdir <- sample(c(-1,0), size=1, prob=probs[-3])
  } else if (Z0 == 0) {
    jmpdir <- sample(c(-1,1), size=1, prob=probs[-2])
  } else {
    jmpdir <- sample(-1:1, size=1, prob=probs)
  }
  jmpposs <- (if (jmpdir == 1) which(!indic.grp)
              else if (jmpdir == 0) which(Zh)
              else which(indic.grp))
  # sample() has a "feature" that it chooses from 1:indic.grp if
  # length(indic.grp)==1.
  if (length(jmpposs) == 1) {
    ind <- jmpposs
  } else {
    ind <- sample(jmpposs, 1)
  }
  if (jmpdir != 0) {
    indic.grp[ind] <- !indic.grp[ind]
    rfwd <- mhIndicGrpJumpMass(ngrp, sgrp, Z0, jmpdir, probs)
    Z0new <- Z0 + (if (is.null(groups) || groups$sizes[ind] > 1) jmpdir else 0)
    rrev <- mhIndicGrpJumpMass(ngrp, sgrp+jmpdir, Z0new, -1*jmpdir, probs)
  }
  list(proposal=indic.grp, dir=jmpdir, ind=ind, rfwd=rfwd, rrev=rrev)
}

# Mass density function for jump proposals for indic.grp
mhIndicGrpJumpMass <- function(ngrp, sgrp, salt, dir, probs) {
  if (dir == 0) {
    NA
  } else if (sgrp == 0) {
    1/ngrp
  } else if (sgrp == ngrp) {
    probs[1]/((1-probs[3])*ngrp)
  } else if (dir == 1 && salt != 0) {
    probs[3]/(ngrp-sgrp)
  } else if (dir == 1) {
    probs[3]/((1-probs[2]) * (ngrp-sgrp))
  } else if (dir == -1 && salt != 0) {
    probs[1]/sgrp
  } else if (dir == -1) {
    probs[1]/((1-probs[2]) * sgrp)
  }
}

updateStateForIndicators <- function(data, state, grpprop, varprop, k) {
  state$indic.grp <- grpprop$proposal
  state$indic.var <- varprop$proposal
  # Only loop over varprop$ind when using multiple traits with indic.grp
  if (is.null(state$indic.grp) || ncol(state$indic.var)==1) {
    state <- updateAuxVars(data, state, k)
  } else {
    for (kk in varprop$ind) {
      state <- updateAuxVars(data, state, kk)
    }
  }
  state
}

# Compute r for proposed indic.grp and indic.var
# Actually compute log(r) until the very end.
rIndicators <- function(params, n, stold, stnew, grpprop, varprop, k) {
  if (is.null(grpprop) || grpprop$dir==0) {
    rln <- log(varprop$rrev) - log(varprop$rfwd)
  } else {
    rln <- log(grpprop$rrev) - log(grpprop$rfwd)
  }
  rln <- rln + computeRIndicatorPrior(params, stold, grpprop, varprop)
  # Only loop over varprop$ind for when using multiple traits with indic.grp
  if (is.null(grpprop) || ncol(stnew$indic.var)==1) {
    rln <- rln + computeRAuxVars(params, n, stold, stnew, k)
  } else {
    for (kk in varprop$ind) {
      rln <- rln + computeRAuxVars(params, n, stold, stnew, kk)
    }
  }
  if (is.null(stnew$OmegaInvDet)) {
    z <- log(1 + n * stnew$tau^2)
    rln <- rln - 0.5 * varprop$dir * length(varprop$ind) * z
  } else {
    rln <- rln - varprop$dir * length(varprop$ind) * log(stnew$tau)
  }
  exp(rln)
}

# Compute contribution to r from the prior for the indicator variables.
computeRIndicatorPrior <- function(params, state, grpprop, varprop) {
  if (params$isOmegaFixed) {
    if (is.null(grpprop)) {
      rln <- sum(log(varprop$omega)) - sum(log(1-varprop$omega))
      rln <- varprop$dir * rln
    } else {
      rln <- log(params$omega.grp) - log(1-params$omega.grp)
      rln <- grpprop$dir * rln
    }
  } else {
    q <- ncol(state$indic.var)
    if (!is.null(grpprop) && grpprop$dir!=0) {
      V <- length(state$indic.grp)
      s <- sum(state$indic.grp)
      dir <- grpprop$dir
      prior.params <- params$omega.grp
    } else if (!is.null(grpprop)) {
      # grpprop$dir == 0
      if (!is.null(params$groups)) {
        inds <- params$groups$group2var[[grpprop$ind]]
        V <- length(inds)
        s <- sum(state$indic.var[inds])
      } else {
        V <- q
        s <- sum(state$indic.var[grpprop$ind,])
      }
      dir <- varprop$dir
      prior.params <- params$omega[grpprop$ind,]
    } else {
      if (q == 1) {
        # basic
        V <- nrow(state$indic.var)
        s <- sum(state$indic.var)
        prior.params <- params$omega
      } else {
        # Multiple traits but not indic.grp
        V <- q
        s <- sum(state$indic.var[varprop$ind,])
        prior.params <- params$omega[varprop$ind,]
      }
      dir <- varprop$dir
    }
    rln <- log(indicatorPriorRatio(V, s, dir, drop(prior.params)))
  }
  rln
}

# Ratio of old and new beta functions reduce to fraction that keeps popping up
# over and over.
indicatorPriorRatio <- function(V, s, dir, params) {
  if (dir == 1) {
    z <- (params["A"] + s) / (params["B"] + V - 1 - s)
  } else if (dir == -1) {
    z <- (params["B"] + V - s) / (params["A"] - 1 + s)
  }
  z
}

### FUNCTIONS FOR tau
updateStateForTau <- function(data, state, tauprop) {
  state$tau <- tauprop
  for (k in seq_along(state$S2)) {
    state <- updateAuxVars(data, state, k)
  }
  state
}

rTau <- function(params, n, stold, stnew) {
  # Ratio of proposal distributions
  rln <- dnorm(stold$tau, mean=stnew$tau, sd=params$tau.sd, log=TRUE)
  rln <- rln - log(pnorm(params$tau.max, mean=stnew$tau, sd=params$tau.sd)
                      - pnorm(params$tau.min, mean=stnew$tau, sd=params$tau.sd))
  rln <- rln - dnorm(stnew$tau, mean=stold$tau, sd=params$tau.sd, log=TRUE)
  rln <- rln + log(pnorm(params$tau.max, mean=stold$tau, sd=params$tau.sd)
                     - pnorm(params$tau.min, mean=stold$tau, sd=params$tau.sd))
  # Ratio of posterior
  for (k in seq_along(stnew$S2)) {
    rln <- rln + computeRAuxVars(params, n, stold, stnew, k)
  }
  svar <- sum(stnew$indic.var)
  if (is.null(stnew$OmegaInvDet)) {
    z <- log(1 + n * stold$tau^2)
    z <- z - log(1 + n * stnew$tau^2)
    rln <- rln + 0.5 * svar * z
  } else {
    rln <- rln + svar * (log(stold$tau) - log(stnew$tau))
  }
  exp(rln)
}
