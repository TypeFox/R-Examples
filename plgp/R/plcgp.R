#******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## lpredprob.CGP:
##
## for the PL resample step -- evaluate the (log) predictive
## density of the new z=(x,c|y) point given suff stats
## and params Zt up to time t by averaging over y following
## Neal's soft-max approach

lpredprob.CGP <- function(z, Zt, prior)
  {
    p <- pred.CGP(z$x, Zt, prior, mcreps=1000, cs=z$c)
    if(!is.finite(p)) warning("bad weight in CGP")    
    return(log(p))
  }


## propagate.CGP:
##
## for the PL propagate step -- add z to Zt and
## calculate the relevant updates to the sufficient
## statistics

propagate.CGP <- function(z, Zt, prior)
  {
    for(i in 2:length(Zt)) {
      
      ## increment t, with sanity check inside pred.GP
      Zt[[i]]$t <- Zt[[i]]$t + 1
      
      ## extract Y
      Yi <- Zt[[i]]$Y; Zt[[i]]$Y <- NULL
    
      ## extend the Y vector in the particle
      tp <- pred.GP(z$x, Zt[[i]], prior, Yi, sub=1:(Zt[[i]]$t-1))
      y <- rt(1, df=tp$df)*sqrt(tp$s2) + tp$m
      if(abs(y) > 100) y <- sign(y)*100
      if(!is.finite(exp(-y))) stop("bad y draw")
      Yi <- c(Yi, y)
      
      ## call the update function on new data, and put new Y
      Zt[[i]]$Y <- Yi      
    }

    ## propose changes to the CGPs
    Zt <- draw.CGP(Zt, prior, l=3, h=4, thin=1)
    
    ## return the propagated particle
    return(Zt)
  }


## prior.CGP:
##
## default prior specification for the CGP model

prior.CGP <- function(m, cov=c("isotropic", "separable", "sim"))
  {
    cov <- match.arg(cov)
    prior <- list(bZero=TRUE, s2p=c(5,40), grate=20, cov=cov)
    if(cov == "isotropic") prior$drate <- 5
    else {
      if(m == 1) stop("use isotropic when m=1")
      prior$drate <- rep(5, m)
      if(cov == "sim") {
        prior$bZero <- TRUE
        prior$drate <- sqrt(1/prior$drate)
      }
    }
    return(prior)
  }


## cv.folds:
##
## calculate a random folds-length partition of indices 1:n

cv.folds <- function (n, folds = 10)
{
  if(n < folds) return(as.list(1:n))
  else return(split(sample(1:n), rep(1:folds, length = n)))
}


## draw.CGP
##
## MH-style draw for the range (d) and nugget (g) parameters
## do the correlation function (K), and the latent (Y) variables

draw.CGP <- function(Zt, prior, l=3, h=4, thin=10)
  {
    ## check if init instead
    if(is.null(Zt)) return(init.CGP(prior))
    
    ## create a Y matrix
    Y <- matrix(NA, nrow=nrow(PL.env$pall$X), ncol=length(Zt))
    Y[,1] <- 0
    
    ## extract each Y & draw parameters of the each GP
    for(i in 2:length(Zt)) { Y[,i] <- Zt[[i]]$Y; Zt$Y <- NULL }

    ## now update latent variables
    for(k in 1:thin) {
      for(j in 2:ncol(Y)) {  ## for each class except the first

        ## loop over folds for block-sampling
        folds <- cv.folds(nrow(Y))
        for(i in 1:length(folds)) {

          ## get the i-th block
          fout <- folds[[i]]; fin <- (1:nrow(Y))[-fout]

          ## calculate the predictive of the new latsents
          tp <- pred.GP(PL.env$pall$X[fout,], Zt[[j]], prior, Y[,j], Sigma=TRUE, sub=fin)
          
          ## calculate the probability of fout classes under the old Ys
          Yj.old <- Y[fout,j]; Yc <- Y[fout,PL.env$pall$C[fout]]
          if(length(fout) > 1) { Yc <- diag(Yc); Ya <- Y[fout,] }
          else Ya <- matrix(Y[fout,], nrow=length(fout))
          pold <- prod(exp(-Yc) / apply(exp(-Ya), 1, sum))

          ## propose new latents
          Yj.prop <- drop(rmvt(1, tp$Sigma, tp$df) + tp$m)
          Yj.prop[abs(Yj.prop) > 100] <- sign(Yj.prop[abs(Yj.prop) > 100])*100
          if(any(!is.finite(exp(-Yj.prop)))) stop("bad Y draw")

          ## calculate probability of fout classes under new Ys
          Y[fout,j] <- Yj.prop; Yc <- Y[fout,PL.env$pall$C[fout]]
          if(length(fout) > 1) { Yc <- diag(Yc); Ya <- Y[fout,] }
          else Ya <- matrix(Y[fout,], nrow=length(fout))
          pnew <- prod(exp(-Yc) / apply(exp(-Ya), 1, sum))

          ## accept or reject accordning to MH
          lalpha <- log(pnew) - log(pold)
          if(runif(1) > exp(lalpha)) Y[fout,j] <- Yj.old
        }
      }
    }

    ## put the latents back in the particles and draw
    for(i in 2:ncol(Y)) {
      ## the extra Y argument signals an updat.GP first
      Zt[[i]] <- draw.GP(Zt[[i]], prior, l=l, h=h, thin=thin, Y=Y[,i])
      Zt[[i]]$Y <- Y[,i]
    }

    ## return the propagated particle
    return(Zt)
  }


## init.CGP:
##
## create a new particle for data X & C, which
## comprises the kitchen-sink (suff stats) variable

init.CGP <- function(prior, d=NULL, g=NULL)
  {
    ## check to make sure all classes are represented
    ## with a coding starting from one
    numclass <- length(unique(PL.env$pall$C))
    if(any(PL.env$pall$C <= 0 || PL.env$pall$C > numclass))
      stop("PL.env$pall$C should range from 1 to numclass")

    ## default d ang if not specified
    if(is.null(d)) d <- 1/prior$drate
    if(is.null(g)) g <- 1/prior$grate
    
    ## initialize numclass-1 regression GPs with
    ## latent variables
    Zt <- list()
    for(i in 2:numclass) {

      ## initial GP parameters
      Zt[[i]] <- list(d=d, g=g)
      Zt[[i]]$t <- nrow(PL.env$pall$X)

      ## initial latent variables
      Y <- 10*(PL.env$pall$C == i)
      Y[Y==0] <- rnorm(sum(Y==0), mean=-10, sd=1)

      ## calculate the sufficient statistics
      Zt[[i]] <- updat.GP(Zt[[i]], prior, Y)
      Zt[[i]]$Y <- Y
    }
      
    ## return the newly allocated particle
    return(Zt)
  }


## pred.CGP:
##
## return the probability of each class at each XX location
## conditional on the Zt particle; the latent y-values
## at the XX locations are integrated out with Monte Carlo
## (mcreps)
##
## this function is optimized for a small XX and large mcreps
## and would need to be re-engineered for large XX and small
## (or one) mcreps

pred.CGP <- function(XX, Zt, prior, mcreps=100, cs=NULL)
  {
    ## coerse the XX input
    XX <- matrix(XX, ncol=ncol(PL.env$pall$X))
    
    ## allocate space for multinomial distn for each XX[i,]
    I <- nrow(XX)
    nclass <- length(Zt)
    if(is.null(cs)) { ## calculating probs for all classes
      p <- data.frame(matrix(NA, nrow=I, ncol=nclass))
      names(p) <- paste("class.", 1:length(Zt), sep="")
    } else p <- rep(NA, I)

    ## set up the matrix of exp(-Y)
    P <- eY <- matrix(NA, nrow=mcreps, ncol=nclass)
    eY[,1] <- 1

    ## store the predictive distribution for every particle
    tp.all <- list()
    for(j in 2:nclass)
      tp.all[[j]] <- pred.GP(XX, Zt[[j]], prior, Zt[[j]]$Y)
    
    ## for each x-location in XX
    for(i in 1:I) {

      ## for each GP (except the first)
      for(j in 2:nclass) {
        
        ## get the parameters to the relevant predictive equations
        tp <- tp.all[[j]][i,]
        
        ## calculate each probability in the multinomal likelihood
        eY[,j] <- exp(-(rt(mcreps, df=tp$df)*sqrt(tp$s2) + tp$m))
        ## eY[!is.finite(eY[,j]),j] <- NA
      }

      ## calculate the probability of each class
      if(!is.null(cs)) p[i] <- mean(eY[,cs]/apply(eY, 1, sum), na.rm=TRUE)
      else {
        ## eY becomes P through softmax, a matrix of class probabilities
        p[i,] <- apply(eY/apply(eY,1,sum,na.rm=TRUE),2,mean,na.rm=TRUE)
      }
    }

    ## return the class probabilities
    return(p)
  }


## params.CGP:
##
## extracts the params from each GP

params.CGP <- function()
  {
    ## extract dimensions
    P <- length(PL.env$peach)
    numGP <- length(PL.env$peach[[1]])-1

    ## allocate data frame (DF) to hold parameters
    params <- data.frame(matrix(NA, nrow=P, ncol=3*numGP))

    ## get the names of the parameters, and set them in the DF
    nam <- c()
    for(i in 2:length(PL.env$peach[[1]])) {
      nam <- c(nam, paste(c("d.", "g.", "lpost."), i, sep=""))
    }
    names(params) <- nam

    ## collect the parameters from the particles
    for(p in 1:P) {
      for(i in 1:numGP) {
        params[p,(i-1)*3+1] <- mean(PL.env$peach[[p]][[i+1]]$d)
        params[p,(i-1)*3+2] <- PL.env$peach[[p]][[i+1]]$g
        params[p,(i-1)*3+3] <- PL.env$peach[[p]][[i+1]]$lpost
      }
    }

    ## return the particles
    return(params)
  }


## latents.CGP:
##
## extracts the latents from each GP

latents.CGP <- function()
  {
    ## not sure yet what to put here.
    P <- length(PL.env$peach)
    Y <- list()
    nam <- paste("Y.", 1:ncol(Y), sep="")
    for(j in 2:length(PL.env$peach[[1]])) {
      Y[[j]] <- data.frame(matrix(NA, nrow=P, ncol=nrow(PL.env$pall$X)))
      names(Y[[j]]) <- nam
    }

    ## collect the latents from the particles
    for(p in 1:P)
      for(j in 2:length(PL.env$peach[[1]]))
        Y[[j]][p,] <- PL.env$peach[[p]][[j]]$Y

    ## return the latents
    return(Y)
  }


## data.CGP:
##
## extract the appropriate columns from the X matrix
## and C vector -- designed to be generic for other cases
## where we would want to get the next observation (end=NULL)
## or a range of observations from begin to end

data.CGP <- function(begin, end=NULL, X, C)
  {
    if(is.null(end) || begin == end)
      return(list(x=X[begin,], c=C[begin]))
    else if(begin > end) stop("must have begin <= end")
    else return(list(x=as.matrix(X[begin:end,]), c=C[begin:end]))
  }


## addpall.CGP:
##
## add data to the pall data structure used as utility
## by all particles

addpall.CGP <- function(Z)
  {
    PL.env$pall$X <- rbind(PL.env$pall$X, Z$x)
    PL.env$pall$C <- c(PL.env$pall$C, Z$c)
    PL.env$pall$D <- distance(PL.env$pall$X)
    if(!is.null(Z$y)) PL.env$pall$Y <- c(PL.env$pall$Y, Z$y)
  }


## entropy.adapt:
##
## return the index into Xcand that has the most potential to
## reduce a high predictive entropy

entropy.adapt <- function(Xcand, rect, prior, verb)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(PL.env$pall$X)+1, " by ME\n", sep="")
    
    ## get predictive distribution information
    outp <- papply(XX=rectscale(Xcand, rect), fun=pred.CGP, prior=prior)

    ## gather the entropy info
    ent <- rep(0, nrow(as.matrix(Xcand)))
    for(i in 1:length(outp)) {
      ent <- ent + calc.ents(outp[[i]])
      ## ent <- ent + drop(apply(outp[[i]], 1, entropy))
      ## ent <- ent + drop(apply(outp[[i]], 1, entropy.bvsb))
    }
    
    ## normalize and return
    return(ent/length(outp))
  }


## data.CGP.adapt:
##
## use the current state of the particules to calculate
## the next adaptive sample from the posterior predictive
## distribution based on the class-entropy

data.CGP.adapt <- function(begin, end=NULL, f, rect, prior, cands=40, verb=2, interp=interp.loess)
  {
    if(!is.null(end) && begin > end) stop("must have begin <= end")
    else if(is.null(end) || begin == end) { ## adaptive sample

      ## choose some adaptive sampling candidates
      if(class(cands) == "function") Xc <- cands()
      else if(is.na(cands)) Xc <- PL.env$Xcand
      else Xc <- lhs(cands, rect)
      ## Xc <- dopt.gp(n=cands, X=NULL, Xc=lhs(10*cands, rect))$XX

      ## calculate the index with the best entropy reduction potential
      as <- entropy.adapt(Xc, rect, prior, verb)
      indx <- which.max(as)
      
      ## return the new adaptive sample
      x <- matrix(Xc[indx,], nrow=1)
      xs <- rectscale(x, rect)

      ## possibly remove the candidate from a fixed set
      if(class(cands) != "function" && is.na(cands)) PL.env$Xcand <- PL.env$Xcand[-indx,]
      
      ## maybe plot something
      if(verb > 1) {
        par(mfrow=c(1,1))
        image(interp(Xc[,1], Xc[,2], as))
        points(rectunscale(PL.env$pall$X, rect))
        points(x[,1], x[,2], pch=18, col="green")
      }

      ## return the adaptively chosen location
      fx <- f(x)
      if(is.list(fx)) return(c(list(x=xs), fx))
      else return(list(x=xs, c=f(x)))

    } else {  ## create an initial design

      ## calculate a LHS 
      ## if(verb > 0) cat("initializing with size", end-begin+1, "LHS\n")
      if(verb > 0) cat("initializing with size", end-begin+1, "MES\n")
      ## X <- lhs(end-begin+1, rect)
      if(class(cands) != "function" && is.na(cands)) Xc <- PL.env$Xcand
      else Xc <- lhs(10*(end-begin+1), rect)
      out <- dopt.gp(end-begin+1, X=NULL, Xcand=Xc)
      X <- out$XX

       ## possibly remove the candidate from a fixed set
      if(class(cands) != "function" && is.na(cands)) PL.env$Xcand <- PL.env$Xcand[-out$fi,]
      
      ## get the class labels
      fX <- f(X)
      if(is.list(fX)) return(c(list(x=rectscale(X, rect)), fX))
      else return(list(x=rectscale(X, rect), c=fX))
    }
  }



## getmap.CGP:
##
## return MAP particle

getmap.CGP <- function(cl=2)
  {
    ## calculate the MAP particle
    mi <- 1
    if(length(PL.env$peach) > 1) {
      for(p in 2:length(PL.env$peach))
        if(PL.env$peach[[p]][[cl]]$lpost > PL.env$peach[[mi]][[cl]]$lpost) mi <- p
    }
    return(PL.env$peach[[mi]])
  }
