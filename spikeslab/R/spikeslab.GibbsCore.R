####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

spikeslab.GibbsCore <- function(n.iter1, n.iter2, 
                         orthogonal, prior,
                         fast, beta.blocks,
                         X, Y, XX, sum.xy,
                         seed,
                         verbose=FALSE,
                         bigp.smalln=FALSE,
                         turn.sigma.on=FALSE,
                         correlation.filter=FALSE,
                         r.eff=NULL,
                         sf= NULL,
                         Wts=NULL,
                         ...)
{
  
  ### --------------------------------------------------------------
  ###
  ###  Global Parameters
  ###  Do not touch unless you know what you are doing !!!
  ###  
  ### --------------------------------------------------------------

  W.shp <- 5                    #cb,mcb
  W.scl <- 50                   #cb,mcb
  lambda.shp <- 2               #lasso:need this to be >=2 for a finite mean
  lambda.scl <- 1               #lasso
  s.shp <- s.scl <- 1e-4        #sigma
  qr.tol <- .Machine$double.eps #qr tolerance
  hv.tolerance <- 1e-5          #used in hypervariance thresholding
  lambda.tolerance <- 1e-10     #used in lasso lambda update
  zap.hypervariance <- 2        #enforcity threshold
  zap.hypervariance.harsh <- 10 #strong enforcity threshold
  factor.bigp.smalln <- 2       #big p small n factor for deciding beta.block
  n.perm <- 10                  #no. permutations for null correlation 
  sparse.tolerance <- 1e-3      #complexity tolerance for numerical stability
  
  ### --------------------------------------------------------------
  ### internal functions

  
  beta.update <- function(b, XX, sum.xy, hyperv, sigma, fast) {
    #------------------------------------------------
    # core draw for beta
    # implements blocked gibbs sampling or
    # straight multivariate normal draw
    #------------------------------------------------
    beta <- b
    nvar <- length(beta)
    if (fast & (nvar >= beta.blocks) ) {
      beta.groups <- beta.fold(nvar, beta.blocks)
      for (j in 1:beta.blocks){
        group.temp <- beta.groups[[j]]
        group.size <- length(group.temp)
        if (length(group.temp) > 1) { 
     	   sum.xy.temp <- sum.xy[group.temp] -
                     as.matrix(XX[group.temp,-group.temp])%*%beta[-group.temp]
         }
         else {
               sum.xy.temp <- sum.xy[group.temp] -
                      t(as.matrix(XX[group.temp,-group.temp]))%*%beta[-group.temp]
         }
	 beta.c.var <- as.matrix(XX[group.temp,group.temp])/sigma +
                           diag(1/(hyperv[group.temp]), nrow=group.size) 
	 qr.var <- qr(beta.c.var, tol = qr.tol)
	 chol.var <- t(chol(beta.c.var))
	 beta[group.temp] <- c(qr.coef(qr.var,
	 chol.var%*%matrix(rnorm(group.size), nrow=group.size) + sum.xy.temp/sigma))
      }
    }
    else {
        beta.c.var <- XX/sigma + diag(1/hyperv, nvar)
	qr.var <- qr(beta.c.var, tol = qr.tol)
        chol.var <- t(chol(beta.c.var))
	beta <- c(qr.coef(qr.var,
	              chol.var%*%matrix(rnorm(nvar), nrow=nvar) + sum.xy/sigma))
    }
    return(beta)
  }
  
  svd.trick <- function(A, G, y, sigma, W, sim = TRUE, beta = NULL) {
    #------------------------------------------------
    # SVD trick for p >> n;  O(p n^2) operation
    # draws beta
    # calculates just GRR (sim=FALSE)
    #------------------------------------------------
    # A (n x p)
    # A = u %*% d %*% t(v)
    # R = u %*% d
    # G (p x p)
    # R (n x n)
    # u (n x n)
    # d (n x n)
    # v (p x n)
    #
    # y (n x 1)     response
    # sigma (p x 1) hypervariances
    # W (p x 1)     complexity.extra
    #
    #
    # ORDER MATTERS FOR COMPUTATIONAL EFFICIENCY !!!
    # no need to update variables with zero cor.xy value
    #
    #
    #
    # Minimum eigenvalue tolerance !!!
    min.eigenv.tol <- 1e-5
    #
    #
    #    
    A <- as.matrix(A)
    n <- nrow(A)
    p <- ncol(A)
    which.nozap <- which(W != 0)
    if (length(which.nozap) == 0) which.nozap <- resample(1:p, 1)
    A <- as.matrix(A[, which.nozap])
    G <- G[which.nozap]
    p.temp <- ncol(A)
    if ((p.temp < n * factor.bigp.smalln) & sim) {#use blocked sampling
      XX.temp <- as.matrix(t(A) %*% A)
      sum.xy.temp <- as.matrix(t(A) %*% y)
      b.return <- rep(0, p)
      b.return[which.nozap] <- beta.update(beta[which.nozap], XX.temp, sum.xy.temp, G, sigma, fast = TRUE)
      return(b.return)   
    }
    else {
      # svd trick: proceed with caution when dealing with small eigenvalues
      # do !!NOT!! use matrix multiplication unless essential (see A.new as an example)
      post.multiply <- pre.multiply <- sqrt(G)
      A.new <- as.matrix(t(t(A) * pre.multiply))
      A.svd <- svd(as.matrix(A.new))
      D <- A.svd$d
#      D[D < min.eigenv.tol] <- min.eigenv.tol
#      D.inv.l <- 1/(D^2+sigma)
#      R <- t(t(A.svd$u)*D)
#      RY <- t(R) %*% y
      D.inv <- D/(D^2 + sigma)
      S.inv <- 1/(D^2 + sigma)
      D.inv[D < min.eigenv.tol] <- 0
      S.inv[D < min.eigenv.tol] <- 0
      V <- A.svd$v
      RU <- t(A.svd$u) %*% y
      beta.c.mean <- V %*% (D.inv * RU)
      if (sim) Z.transform <- sqrt(S.inv) * rnorm(n) 
      if (correlation.filter) {
        b.return <- rep(0, p)
        if (sim) {
          b.return[which.nozap] <- post.multiply * (beta.c.mean + sqrt(sigma) * (V %*% Z.transform))
        }
        else {
          b.return[which.nozap] <- post.multiply * (beta.c.mean)
        }
        b.return
      }
      else {
        if (sim) {
          post.multiply * (beta.c.mean + sqrt(sigma) * (V %*% Z.transform))
        }
        else {
          post.multiply * (beta.c.mean)
        }
      }
    }
  }



  ### --------------------------------------------------------------
  ###
  ### initialize Gibbs parameters (same for all priors)
  ### calibrate V.small to be consistent across priors
  ###
  ### --------------------------------------------------------------

  #dimensions; vectors; arrays
  n.cov <- ncol(X)
  n.data <- length(Y)
  b.m <- b.percent <- n.sample <- 0
  complexity.vec <- lambda.vec <- sigma.vec <- rep(0, n.iter2)
  b.gamma <- model <- vector("list", n.iter2)
  
  #hyperparameters
  if (orthogonal) V.small <- 1e-5 else V.small <- 5e-3
  calibrate.cb <- calibrate.bimodal <- calibrate.mcb <- W.scl/(W.shp-1)
  calibrate.lasso <- (W.scl/(W.shp-1))*(lambda.shp-1)/(2*lambda.scl)

  #initialize parameters
  sigma <- n.data
  beta <- rnorm(n.cov, sd=0.25)
  tau <- rgamma(n.cov, 1)	
  Waugment <- Vaugment <- rho <- rep(1, n.cov)
  complexity <- 0.5
  complexity.extra <- rep(1, n.cov)
  if (correlation.filter) {
    prior.type <- 0
    V.small <- sqrt(V.small)
    if (bigp.smalln & is.null(r.eff)) {
      cor.xy <- c(atanh(sum.xy/sqrt(n.data*sum(Y^2))))
      sum.xy.perm <- 0
      for (b in 1:n.perm) {
          sum.xy.perm <- sum.xy.perm + t(X[resample(1:n.data, n.data, replace = F), ])%*%Y
      }
      sum.xy.perm <- sum.xy.perm/n.perm
      cor.permute.xy <- c(atanh(sum.xy.perm/sqrt(n.data*sum(Y^2))))
      cor.cut.off1 <- 1.5*sd(cor.permute.xy)
      cor.cut.off2 <- 2*sd(cor.permute.xy)
      cor.cut.off3 <- 2.5*sd(cor.permute.xy)
      cor.cut.off4 <- 3*sd(cor.permute.xy)
      pt.0 <- (abs(cor.xy) <= cor.cut.off1)
      pt.1 <- (cor.cut.off1 < abs(cor.xy) & abs(cor.xy) <= cor.cut.off2)
      pt.2 <- (cor.cut.off2 < abs(cor.xy) & abs(cor.xy) <= cor.cut.off3)
      pt.3 <- (cor.cut.off3 < abs(cor.xy) & abs(cor.xy) <= cor.cut.off4)
      pt.4 <- (cor.cut.off4 < abs(cor.xy))
      if (is.null(r.eff)) {
        r.eff <- vector("list", 5)
        r.eff[[1]] <- which(pt.0)
        r.eff[[2]] <- which(pt.1)
        r.eff[[3]] <- which(pt.2)
        r.eff[[4]] <- which(pt.3)
        r.eff[[5]] <- which(pt.4)
      }
    }
    else {
      #under pre-assigned groupings, complexity is multiplied using data-driven values
      #UNCOMMENT THE FOLLOWING FOR THIS EFFECT
      #cor.xy <- c(atanh(sum.xy/sqrt(n.data*sum(Y^2))))
      #cor.permute.xy <- apply(cbind(1:n.cov), 1, function(i){atanh(Cor(Y, resample(X[,i])))})
      #complexity.extra <- apply(cbind(cor.xy), 1, function(r) {mean(abs(cor.permute.xy) < abs(r))})
    }
  }
  else {
    prior.type <- match(prior, c("cb", "mcb", "lasso"))
  }
  hypervariance <- rho*tau

  #prior specific details
  if (prior == "cb") {
    lambda <- 0
  }
  else if (prior == "mcb") {
    #!!!calibration for mcb not correct yet!!!
    #exp.scl <- sqrt(2/(V.small*calibrate.mcb))
    exp.scl <- 3
    V.small <- lambda <- 0
    
  }
  else if (prior == "bimodal") {
    lambda <- 0
    V.small <-  V.small*calibrate.bimodal
    V.big <- n.data
    hypervariance <- resample(c(V.small,V.big), size=n.cov, replace=TRUE)
  }
  else {
    V.small <- V.small*calibrate.lasso 
    lambda <- rgamma(1,2,1)
  }

  # random effects
  # if f.eff NULL  ---> shared complexities + usual shrinkage on r effects
  # if f.eff !NULL ---> shared complexities + ZERO shrinkage on f/r effects (not enforced)
  f.eff <- NULL
  if (!is.null(r.eff)) {
    complexity.reff <- rep(0.5, length(r.eff))
    f.eff <- setdiff(1:n.cov, unlist(r.eff))
    if (length(f.eff) == 0) {
      complexity.eff.vec <- matrix(0, n.iter2, length(r.eff))
      complexity.feff <- f.eff <- NULL
    }
    else {
     complexity.eff.vec <- matrix(0, n.iter2, length(r.eff) + 1)
      complexity.feff <- 0.5
    }
  }

  ### ---------------------------------
  #  Pretty terminal output

  if (verbose) {
    cat("                                                                                    \r")
  }


  ### -------------------------------------------------------------
  ###
  ### Burn in followed by sampled values
  ###
  ### -------------------------------------------------------------

  for (i in 1:(n.iter1+n.iter2)){

        ### ---------------------------------	
        # sample beta
        # remember that sigma is scaled by n

        if (orthogonal) {
           beta.c.sd <- (n.data/sigma + 1/hypervariance)**(-0.5)
           beta.c.m <- (beta.c.sd^2)*sum.xy/sigma
           beta <- rnorm(n.cov, mean=beta.c.m, sd=beta.c.sd)
        }
        else {
          if (bigp.smalln) {#SVD trick for big p small n
            beta <- svd.trick(X, hypervariance, Y, sigma, complexity.extra, beta = beta)
          }
          else {#beta-block update
            beta <- beta.update(beta, XX, sum.xy, hypervariance, sigma, fast)
          }
        }
        
        #save sparse beta
        if (i > n.iter1) {
            if (bigp.smalln) zap.hvar <- zap.hypervariance.harsh else zap.hvar <- zap.hypervariance
            hv.pt <- (hypervariance > zap.hvar)
            hvar <- hypervariance
          if (orthogonal & bigp.smalln) {
            beta.sparse <- beta
            beta.sparse[!hv.pt] <- 0
          }
          else {
            beta.sparse <- rep(0, n.cov)
            if (sum(hv.pt) > 0) {
              if (sum(hv.pt) < n.data & !bigp.smalln) {
                beta.c.var <- as.matrix(XX[hv.pt, hv.pt])/sigma +
                  diag(1/hvar[hv.pt], nrow=sum(hv.pt))
                qr.var <- qr(beta.c.var, tol = qr.tol)
                beta.sparse[hv.pt] <- qr.coef(qr.var, sum.xy[hv.pt]/sigma)
              }
              else {
                if (sum(hv.pt) < n.data & bigp.smalln) {
                  XX.hv <- t(as.matrix(X[ ,hv.pt]))%*%(as.matrix(X[ ,hv.pt]))
                  beta.c.var <- as.matrix(XX.hv/sigma + diag(1/hvar[hv.pt], nrow=sum(hv.pt)))
                  qr.var <- qr(beta.c.var, tol = qr.tol)
                  beta.sparse[hv.pt] <- qr.coef(qr.var, sum.xy[hv.pt]/sigma)
                }
                else {#use SVD trick
                  beta.sparse[hv.pt] <- svd.trick(X[, hv.pt], hvar[hv.pt], Y,
                                        sigma, complexity.extra[hv.pt], sim = FALSE)
                }
              }
            }
          }
        }

        ### ---------------------------------	
        # sample hypervariance
        # all fixed effects, or fixed + random effects
        
        if (is.null(r.eff)) {

          ### FIXED EFFECTS
          
          #...sample gamma under N(0, W{-1}Gamma) prior
          #suffices to scale betas by sqrt(W_p)
          if (!is.null(Wts)) {
            beta.temp <- beta
            beta <- sqrt(Wts) * beta
          }

          if (prior != "bimodal") {
            #if (prior == "lasso") {
            #  tau <- 1/rinvgauss(n.cov , lambda/(pmax(sqrt(rho)*abs(beta) , lambda.tolerance)),
            #                   lambda^2/rho) 
            #}
            #if (prior == "mcb") {
            #   Waugment <- 1/rinvgauss(n.cov , exp.scl/(pmax(abs(beta) , lambda.tolerance)), exp.scl^2)
            #   Waugment[rho != V.small] <- rexp(sum(rho != V.small), rate=exp.scl^2/2)
            #}
            spikeSlabVar <- .C("spikeSlabVar",
                          as.double(beta),
                          W=as.double(tau),
                          as.double(Waugment),
                          V=as.double(rho),
                          V.extra=as.double(Vaugment),
                          as.integer(prior.type),
                          as.integer(n.cov),
                          as.double(W.scl),
                          as.double(W.shp),
                          as.double(V.small),
                          as.double(complexity),
                          as.double(complexity.extra),
                          seed=as.integer(seed))
            tau <- spikeSlabVar$W
            rho <- spikeSlabVar$V
            seed <- spikeSlabVar$seed
            if (correlation.filter) {#cor.xy filter details
              rho.extra <- spikeSlabVar$V.extra
              hypervariance <- tau*rho*rho.extra
              hypervariance[which(complexity.extra == 0)] <- hv.tolerance#cor.xy=0  
            }
            else {
              if (prior != "mcb") {
                hypervariance <- tau*rho
              }
              else {
                hypervariance <- tau
              }
            }
          }
          else {#two-component prior
             rho <- hypervariance <- apply(cbind(beta),1,varianceSample,w=complexity,v0=V.small,V=V.big)
          }

          #...now scale beta back if using a N(0, W^{-1}Gamma) prior
          if (!is.null(Wts)) {
            beta <- beta.temp
          }

        }
        else {###FIXED + RANDOM EFFECTS
          
          ### RANDOM EFFECTS
          for (k in 1:length(r.eff)) {
            r.eff.k <- r.eff[[k]]
            if (!all(complexity.extra[r.eff.k] == 0)) {
              spikeSlabVar <- .C("spikeSlabVar",
                          as.double(beta[r.eff.k]),
                          W=as.double(tau[r.eff.k]),
                          as.double(Waugment[r.eff.k]),
                          V=as.double(rho[r.eff.k]),
                          V.extra=as.double(Vaugment[r.eff.k]),
                          as.integer(prior.type),
                          as.integer(length(r.eff.k)),
                          as.double(W.scl),
                          as.double(W.shp),
                          as.double(V.small),
                          as.double(complexity.reff[k]),
                          as.double(complexity.extra[r.eff.k]),
                          seed=as.integer(seed))
              tau[r.eff.k] <- spikeSlabVar$W
              rho[r.eff.k] <- spikeSlabVar$V
              seed <- spikeSlabVar$seed
              if (correlation.filter) {#cor.xy filter details
                hypervariance[r.eff.k] <- tau[r.eff.k]*rho[r.eff.k]*spikeSlabVar$V.extra
              }
              else {
                hypervariance[r.eff.k] <- tau[r.eff.k]*rho[r.eff.k]
              }
            }
          }
          
          ### FIXED EFFECTS
          if (!is.null(f.eff)) {
            spikeSlabVar <- .C("spikeSlabVar",
                          as.double(beta[f.eff]),
                          W=as.double(tau[f.eff]),
                          as.double(Waugment[f.eff]),
                          V=as.double(rho[f.eff]),
                          V.extra=as.double(Vaugment[f.eff]),
                          as.integer(prior.type),
                          as.integer(length(f.eff)),
                          as.double(W.scl),
                          as.double(W.shp),
                          as.double(V.small),
                          as.double(complexity.feff),
                          as.double(complexity.extra[f.eff]),
                          seed=as.integer(seed))
            tau[f.eff] <- spikeSlabVar$W
            rho[f.eff] <- spikeSlabVar$V
            seed <- spikeSlabVar$seed
            if (correlation.filter) {#cor.xy filter details
                hypervariance[f.eff] <- tau[f.eff]*rho[f.eff]*spikeSlabVar$V.extra
            }
            else {
                hypervariance[f.eff] <- tau[f.eff]*rho[f.eff]
            }
          }
          rho[which(complexity.extra == 0)] <- V.small#cor.xy=0          
          hypervariance[which(complexity.extra == 0)] <- hv.tolerance#cor.xy=0
        }

        ### ---------------------------------	
        # sample complexity

        if (is.null(r.eff)) {
          complexity <- rbeta(1,(1 + sum(1*(rho != V.small))), (1 + sum(1*(rho == V.small))))
          if (1 - complexity < sparse.tolerance) complexity <- 1 - sparse.tolerance
          if (complexity < sparse.tolerance) complexity <- sparse.tolerance
        }
        else {
          for (k in 1:length(r.eff)) {
            rho.k <- rho[r.eff[[k]]]
            complexity.k <- complexity.reff[k] <- rbeta(1,(1 + sum(1*(rho.k != V.small))),
                                                        (1 + sum(1*(rho.k == V.small))))
            if (1 - complexity.k < sparse.tolerance) complexity.k <- 1 - sparse.tolerance
            if (complexity.k < sparse.tolerance) complexity.k <- sparse.tolerance
            complexity.reff[k] <- complexity.k
          }
          if (!is.null(f.eff)) {
            rho.k <- rho[f.eff]
            complexity.feff <- rbeta(1,(1+sum(1*(rho.k != V.small))),
                                                       (1+sum(1*(rho.k == V.small))))
            if (1 - complexity.feff < sparse.tolerance) complexity.feff <- 1 - sparse.tolerance
            if (complexity.feff < sparse.tolerance) complexity.feff <- sparse.tolerance
          }
        }

        ### ---------------------------------	
        # sample lambda
        
	if (prior == "lasso") {
	  lambda <- sqrt(rgamma(1, lambda.shp+sum(tau)/2)/(lambda.scl+n.cov))
          lambda <- max(1e-6, lambda)
        }

        ### ---------------------------------	
        # sample sigma
        # note that sigma is scaled by n 

        if (turn.sigma.on) {
          if (is.null(sf)) {
             sigma <- n.data * (s.scl+sum((Y-(X%*%beta))^2)/(2*n.data)) /
	                     rgamma(1,s.shp+n.data/2)
          }
          else {
             sigma <- n.data * (s.scl+sum((Y/sf-(X%*%beta))^2)/(2*n.data)) /
	                     rgamma(1,s.shp+n.data/2)
          }
        }
      
        ### ---------------------------------			
        # verbose details
        # save output values
        
        if (i<=n.iter1 & i%%50==0 & verbose) {
          cat("                                                                                    \r")
          if (is.null(r.eff)) {
            cat("\t \t", i, ":", round(complexity, 5), "\r")
          }
          else {
            cat("\t \t", i, ":", round(c(complexity.feff, complexity.reff), 5), "\r")
          }
        }
        if (i>n.iter1){
          n.sample <- n.sample + 1
	  b.m <- b.m + beta.sparse
          b.percent <- b.percent + 1 * (rho != V.small)
	  complexity.vec[n.sample] <- complexity
          if (!bigp.smalln) {
            b.gamma[[n.sample]] <- sigma  / hypervariance
            model[[n.sample]] <- intersect(which(rho != V.small), which(beta.sparse != 0))
          }
          ### f/r effects complexities
          if (!is.null(r.eff)) {
            if (is.null(f.eff)) {
              complexity.eff.vec[n.sample, ] <- complexity.reff
            }
            else {
              complexity.eff.vec[n.sample, ] <- c(complexity.reff,complexity.feff)
            }
          }
          lambda.vec[n.sample] <- lambda
          sigma.vec[n.sample] <- sigma/n.data
          m.size <- sum(rho != V.small)
	  if (i%%50==0 & verbose) {
            cat("                                                                                    \r")
            if (is.null(r.eff)) {
              cat("\t \t", n.sample, ":", round(complexity, 5), "\r")
            }
            else {
              cat("\t \t", n.sample, ":", round(c(complexity.feff, complexity.reff), 5), "\r")
            }
          }
        }
  }

  ### ------------------------------------
  ### decide what complexity to return
  if (!is.null(r.eff)) complexity.vec <- complexity.eff.vec
  ### ------------------------------------
  

  ### -------------------------------------------------------------
  ###
  ### Return the goodies
  ###
  ### -------------------------------------------------------------

  out <- list(
       b.m=b.m/n.sample,
       b.percent=b.percent/n.sample,
       complexity.vec=complexity.vec,
       lambda.vec=lambda.vec,
       sigma.vec=sigma.vec,
       hyperv=b.gamma,
       model=model)

  return(out)
  
}
