#### FIXME: hay que arreglar trace.plot!!!!!




####  Copyright (C) 2005-2008  Oscar Rueda Palacio and Ramon Diaz-Uriarte
#### This program is free software; you can redistribute it and/or
#### modify it under the terms of the GNU General Public License
#### as published by the Free Software Foundation; either version 2
#### of the License, or (at your option) any later version.

#### This program is distributed in the hope that it will be useful,
#### but WITHOUT ANY WARRANTY; without even the implied warranty of
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#### GNU General Public License for more details.

#### You should have received a copy of the GNU General Public License
#### along with this program; if not, write to the Free Software
#### Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
#### USA.


### zz: FIXME: remove .__DELETE_GZIPPED <- TRUE ## turn to true for real
####            .__DELETE_BETA <- FALSE  ## later set to TRUE, when we have summaries from C
####    return name of temporary directory!! and have a function for extraction

#### TODO:
####      - Rewrite getEdges using the stored viterbi sequences!!!!

#### DEFINE SUPERSEED; only needed for debugging
#### SUPERSEED <- 123

#### Includes the code that returns probs of states being gained/lost/no-change


#### For now until we set in stone what we do with dummy chroms.
##   this is just for having nice names in output

dummyChr <- "-"


.__DELETE_GZIPPED <- TRUE ## turn to true for real
.__RJACGH_DEBUG <- FALSE
.__DELETE_BETA <- FALSE  ## later set to TRUE, when we have summaries from C

#########################################################################
##Functions for Non Homogeneous Hidden Markov Models with normal
##probability emissions
#########################################################################

#########################################################################
##Function to compute the transition matrix
## q are the intercept probabilities
## beta is the vector of regression coefficients
## x is the distance to the next gene
#########################################################################
QNH <- function(beta, x, q=-beta) {

  k <- nrow(q)
  Q <- exp(q + beta*x)
  Q <- Q / rowSums(Q)
  Q
}

#########################################################################
## filter and likelihood of the non - homogeneous HMM
## q, beta, x are the same as before
## mu & sigma.2 are the normal distribution parameters
## stat is the initial distribution
## k is the number of hidden states
#########################################################################

normal.HMM.likelihood.NH.filter <- function(y, k, beta, x, mu, sigma.2,
                                     stat=NULL, q=-beta) {
  ## Prevents that c, denominator never gets zero
  TOL <- 1E-10
  if (is.null(stat)) stat <- rep(1/k, k)
  
  ##last x is not observed. Besides, it doesn't count on the likelihood
  x <- c(x, 0)
  
  n <- length(y)
  if (k==1) {
    loglik <- sum(dnorm(y, mu, sqrt(sigma.2), log=TRUE))
    filter.cond <- matrix(rep(1, n+1), n+1)
    filter <- matrix(rep(1, n), n)
    res <- list(loglik=loglik, filter.cond=filter.cond, filter=filter)
  }
  else {
    filter <- matrix(NA, n,k)
    filter.cond <- matrix(NA, n+1,k)
    c <- rep(NA, n)
    filter.cond[1,] <- stat
    for (i in 1:n) {
      c[i] <- sum(filter.cond[i,]*dnorm(y[i], mu, sqrt(sigma.2)))
      c[i] <- ifelse(identical(c[i], 0), TOL, c[i])
      filter[i,] <- filter.cond[i,]*dnorm(y[i], mu, sqrt(sigma.2))/c[i]
      Q <- QNH(q=q, beta=beta, x=x[i])
      for (j in 1:k) {
        filter.cond[i+1,j] <- sum(filter[i,]*Q[,j,drop=FALSE])
      }
    }
    res <- list(loglik=sum(log(c)), filter.cond=filter.cond, filter=filter, c=log(c))
  }
  res
}

#########################################################################
## likelihood of the non - homogeneous HMM. C version
## q, beta, x are the same as before
## mu & sigma.2 are the normal distribution parameters
## stat is the initial distribution
## k is the number of hidden states
#########################################################################

normal.HMM.likelihood.NH.C <- function(y, x, mu,
                                       sigma.2, beta, stat=NULL) {
  k <- length(mu)
  if (is.null(stat)) stat <- rep(1/k, k)
  ##last x is not observed. Besides, it doesn't count on the likelihood
  x <- c(x, 0)
  q <- -beta
  loglik <- .C("normalNHHMMlikelihood", y=as.double(y), k=as.integer(k),
               x=as.double(x), n=as.integer(length(y)),
               q=as.double(as.vector(q)), ## beta=as.double(beta),
               stat=as.double(stat), mu=as.double(mu),
               sigma.2=as.double(sigma.2),
               loglik=as.double(0))
  loglik$x <- loglik$x[-length(x)]
  loglik
  
}


#########################################################################
## Function to simulate hidden state sequence by backward smoothing
## The call to normal.HMM.likelihood.NH should be changed for c++ version
## Slightly modified. Doesn't sample, it chooses the state with maximal
## probability.
## now can do model averaging
## now viterbi is used instead of this one
#########################################################################
hidden.states <- function(obj, x) {
  k <- length(obj$mu)
    n <- length(obj$y)
    filter <- normal.HMM.likelihood.NH.filter(y=obj$y, k=k, q=-obj$beta,
                                       beta=obj$beta, x=x, mu=obj$mu,
                                       sigma.2=obj$sigma.2, stat=obj$stat)$filter
    states <-rep(NA, n)
    B <- array(NA, c(n, k))
    B[n,] <- filter[n,]
    states[n] <- which.max(B[n,])
    for (i in (n-1):1) {
      Q <- QNH(q=-obj$beta, beta=obj$beta, x=x[i])
      den <- sum(filter[i,]*Q[,states[i+1]])
      den <-ifelse(den==0, 1, den)
      B[i,] <- filter[i,]*Q[,states[i+1]]
      B[i,] <- B[i,] / den
      states[i] <- which.max(B[i,])
    }
  res <- list(states=states, prob.states=B)
  res
}

simulateRJaCGH <- function(n, x=NULL, mu, sigma.2, beta, start) {
  y <- rep(NA, n)
  states <- rep(NA, n)
  k <- length(mu)
  if (is.null(x)) x <- rep(0, n-1)
  states[1] <- start
  y[1] <- rnorm(1, mu[states[1]], sqrt(sigma.2[states[1]]))
  for (i in 2:n) {
    Q <- QNH(beta=beta, x=x[i-1], q=-beta)
    states[i] <- sample(x=1:k, size=1, prob=Q[states[i-1],])
    y[i] <- rnorm(1, mu[states[i]], sqrt(sigma.2[states[i]]))
  }
  list(states=states, y=y)
}

#########################################################################
## Plot transition probabilities
#########################################################################

plotQNH <- function(x, beta, q=-beta, col=NULL,...) {
  n<- length(x)
  k <- nrow(q)
  prob <- matrix(NA, n, k)
  if (is.null(col)) {
    col <- rep(1, k)
    col[grep("[G]", rownames(beta))] <- 2
    col[grep("[L]", rownames(beta))] <- 3
  }
  for (i in 1:n) {
    prob[i,] <- diag(QNH(q=q, beta=beta, x=x[i]))
  }
  plot(prob[order(x),1] ~ sort(x), type="l", ylim=c(min(prob), max(prob)), col=col[1], ...)
  ## case with homogeneous probabilities
  if (min(x)==0 && max(x)==0) abline(h=prob[1,1], col=col[1])
  if (k >1) {
    for (i in 2:k) {
      lines(prob[order(x),i] ~ sort(x), col=col[i])
        ## case with homogeneous probabilities
      if (min(x)==0 && max(x)==0) abline(h=prob[1,i], col=col[i])

    }
  }
  rug(x)
}

 
get.viterbi.seq2 <- function(res) {
    ffnames <- unlist(strsplit(res$viterbi.filename, "\n"))
    ## the following does not work well
###     mapply(function(x, y, z)
###            list(num_sequences = x,
###                 sum_mcmc_iter = y,
###                 gzipped_filename = z),
###            x = res$viterbi.num_sequences,
###            y = res$viterbi.sum_mcmc_iter,
###            z = ffnames)

    retv <- list()
    n1 <- length(res$viterbi.num_sequences)
    if(n1 != length(ffnames))
        stop("FATAL error in get.viterbi.seq2")
###     if(n1 == 1)
###         return(list(num_sequences = res$viterbi.num_sequences[1],
###                     sum_mcmc_iter = res$viterbi.sum_mcmc_iter,
###                     gzipped_filename = ffnames[1]))
###     else {
    for(i in 1:n1) {
        if(!file.exists(ffnames[i]))
            cat("\n get.viterbi.seq2 cannot find file ", ffnames[i]," \n")
        retv[[i]] <- list()
        retv[[i]]$num_sequences <- res$viterbi.num_sequences[i]
        retv[[i]]$sum_mcmc_iter <- res$viterbi.sum_mcmc_iter
        retv[[i]]$gzipped_filename <- ffnames[i]
    }
    return(retv)
###     }
}
 

MetropolisSweep.C <- function(y, x, k.max, Chrom, model=NULL,
                              var.equal,
                              maxVar=diff(range(y)),
                              normal.reference, 
                              burnin, TOT, prob.k, pb, ps, mu.alfa,
                              mu.beta, s1, s2, init.mu,
                              init.sigma.2,
                              init.beta, sigma.tau.mu, sigma.tau.sigma.2,
                              sigma.tau.beta, tau.split.mu=NULL,
                              tau.split.beta=NULL, stat, start.k,
                              RJ=TRUE,
                              NC, deltaT,
                              write_seq,
                              window = NULL,
                              singleState = FALSE,
                              delete_gzipped = .__DELETE_GZIPPED,
                              temp.dir) {
  n <- length(y)
  ##Size of vectors
  ## Now we've had to put 2 * TOT beacuse of swap move
  ## In a sweep, we can get two values for a same r
  size.mu <- 2 * TOT * k.max*(k.max+1)/2
  size.sigma.2 <- 2 * TOT * k.max*(k.max+1)/2
  size.beta <- 2 * TOT * k.max * (k.max+1) * (2*k.max+1) / 6
  mu <- rep(0, size.mu)
  sigma.2 <- rep(0, size.sigma.2)
  beta <- rep(0, size.beta)
  probStates <- rep(0, n*(k.max^2 - k.max) / 2)
  loglik <- rep(0, 2 * TOT * k.max)
  if (is.null(start.k)) start.k <- 0

  ## new parameter for split/combine
  tau.split.beta <- 1/ (1-median(x))
  
  if(model != "Chrom") {
      index <- c(which(!duplicated(Chrom)) - 1, length(y))
      genome.param <- length(index) - 1
      num_sequences <- rep(-99, genome.param)
  } else { 
      index <- c(0, length(y))
      genome.param <- 1
      num_sequences <- -99
  }
  if(length(y) != (length(x) + 1)) {
      stop(paste("Length of vector of distances and data vector are different in MetropolisSweep.C", model))
  }
  ## viterbi seqs stuff
  if(write_seq) {
      filename <- tempfile2(pattern = "rjacgh_seq", tmpdir = temp.dir)
      if(model != "Chrom") {
          filenames <- paste(filename, "_chr_", unique(Chrom),
                             ".gz", sep = "")
          lapply(filenames, file.create) ## create now, so no same name later
          filename <- paste(filenames, collapse="\n")
      } else {
          if(is.null(Chrom)) {
              chtmp <- "na"
          } else {
              chtmp <- Chrom[1]
          }
          filename <- paste(filename, "_chr_", chtmp, ".gz", sep = "")
          file.create(filename)
      }
  } else {
      filename <- NULL
  }
## set.seed(SUPERSEED)
## cat("\n            RUNIF BEFORE METROPOLIS ", runif(1), "\n")
  res <- .C("MetropolisSweep", y=as.double(y), x=as.double(x),
            varEqual=as.integer(var.equal), genome= as.integer(genome.param),
            index=as.integer(index),
            kMax=as.integer(k.max), n=as.integer(length(y)),
            burnin=as.integer(burnin),
            TOT=as.integer(TOT), times=as.integer(rep(0, k.max)),
            probB=as.integer(c(0,0)), probD=as.integer(c(0,0)),
            probS=as.integer(0),
            probC=as.integer(0),
            probE=as.integer(c(0,0)),
            probK=as.double(prob.k),
            pb=as.double(pb), ps=as.double(ps),
            muAlfa=as.double(mu.alfa),
            muBeta=as.double(mu.beta),
            s1=as.double(s1),
            s2=as.double(s2), init.mu=as.double(init.mu),
            init.sigma.2=as.double(init.sigma.2),
            init.beta=as.double(init.beta), 
            sigmaTauMu=as.double(sigma.tau.mu),
            sigmaTauSigma2=as.double(sigma.tau.sigma.2),
            sigmaTauBeta=as.double(sigma.tau.beta),
            tauSplitMu=as.double(tau.split.mu),
            tauSplitBeta=as.double(tau.split.beta),
            k=as.integer(rep(0, 3*(TOT-1)+1)), mu=as.double(mu),
            sigma2=as.double(sigma.2), beta=as.double(beta),
            stat=as.double(stat), startK=as.integer(start.k),
            RJ=as.integer(1*RJ), maxVar=as.double(maxVar),
            probStates=as.double(probStates),
            loglik=as.double(loglik),
            NC=as.integer(NC),
            deltaT=as.double(deltaT),
            write_seq = as.integer(write_seq),
            filename = as.character(filename),
            num_sequences = as.integer(num_sequences)
            )
## cat("\n            RUNIF AFTER METROPOLIS ", runif(1), "\n")
  ##Reconstruct objects
  obj <- list()
  obj$fit.k <- list()
  gc()
  indexStat <- 1
  indexStates <- 1
  indexJointStates <- 1
  indexMu <- 1
  indexBeta <- 1
  indexLoglik <- 1
  for (i in 1:k.max) {
    obj$fit.k[[i]] <- list()
    obj$fit.k[[i]]$stat <- res$stat[indexStat:(indexStat + i-1)]
    ## we create some tmps so that computing prob.xyz does not require readRO2
    tmp.mu <- matrix(res$mu[indexMu: (indexMu + res$times[i]*i -1)],
                     ncol=i, byrow=TRUE)
    tmp.sigma.2 <- matrix(res$sigma2[indexMu: (indexMu + res$times[i]*i -1)],
                          ncol=i, byrow=TRUE)
    tmp.loglik <- res$loglik[indexLoglik:(indexLoglik + res$times[i] -1)]
    tmp.beta <- array(res$beta[indexBeta: (indexBeta + res$times[i]*i*i-1)],
                      dim=c(i, i, res$times[i]))

    obj$fit.k[[i]]$mu <- writeRO2(tmp.mu, "mu", temp.dir)
    obj$fit.k[[i]]$sigma.2 <- writeRO2(tmp.sigma.2, "sigma.2", temp.dir)

    obj$fit.k[[i]]$loglik <- writeRO2(tmp.loglik, "loglik", temp.dir)
    obj$fit.k[[i]]$beta <- writeRO2(tmp.beta, "beta", temp.dir)

    indexMu <-  indexMu + TOT*i
    indexLoglik <- indexLoglik + TOT
    indexBeta <- indexBeta + TOT*i*i
    indexStat <- indexStat + i

    frull <- function(x) return(length(unique(x))/length(x))
    if(i == 1) {
        ## zz FIXME: all except beta are the same for i==1 or more!!!
        obj$fit.k[[1]]$prob.mu <- frull(tmp.mu)
        obj$fit.k[[1]]$prob.sigma.2 <- frull(tmp.sigma.2)
        obj$fit.k[[1]]$prob.beta <- frull(tmp.beta) ## but this is always 1! zz??

        obj$fit.k[[i]]$prob.states <- writeRO2(matrix(rep(1, length(y)), ncol=1),
                                               "prob.states", temp.dir)

    } else {
        obj$fit.k[[i]]$prob.mu <- frull(tmp.mu[, 1])
        obj$fit.k[[i]]$prob.sigma.2 <- frull(tmp.sigma.2[, 1])
        obj$fit.k[[i]]$prob.beta <- frull(tmp.beta[2, 1, ])

        if (nrow(tmp.mu) > 0) {
            tmp.psts <- res$probStates[indexStates : (indexStates +
                                                      n*(i-1) -1)]
            tmp.psts <- matrix(tmp.psts, nrow=n, ncol=i-1)
            obj$fit.k[[i]]$prob.states <- writeRO2(cbind(tmp.psts, 1 - rowSums(tmp.psts)),
                                                   "prob.states", temp.dir)
        }
        else {
            obj$fit.k[[i]]$prob.states <- NULL
        }
    }
    indexStates <- indexStates + n*(i-1)
  }

  obj$y <- writeRO2(y, "y", temp.dir)
                    ##it has to be here, as simpler for many calls.
  obj$sdY <- sd(y)
  obj$lengthY <- length(y)
  obj$k <- res$k
  ## If there are random moves we won't have 3*TOT k values ##
  ## we take out the missing values (zero's) ##
  obj$k <- obj$k[obj$k > 0]
  obj$k <- factor(obj$k , levels=1:k.max)
  
  ## Relabel states
  sllist <- relabel.core(obj = obj,
                         normal.reference = normal.reference,
                         window = window,
                         singleState = singleState)
  k.max <- max(as.numeric(levels(obj$k)))
  for(i in (1:k.max))
      obj$fit.k[[i]]$state.labels <- sllist[[i]]
  
  obj$prob.b <- res$probB
  obj$prob.d <- res$probD
  obj$prob.s <- res$probS
  obj$prob.c <- res$probC
  obj$prob.e <- res$probE
##   obj$x <- x
##   obj$x[obj$x==-1] <- NA
  obj$viterbi.filename <- filename
  obj$viterbi.num_sequences <- res$num_sequences
  obj$viterbi.sum_mcmc_iter <- sum(unlist(res$times))
  
  rm(res)
  gc()
  class(obj) <- "internal.Metropolis.Sweep.C"
  obj
}


RJMCMC.NH.HMM.Metropolis <- function(y, Chrom=NULL, x=NULL,
                                     index=NULL, maxVar, 
                                     model=NULL, var.equal, max.dist=NULL,
                                     normal.reference=0, 
                                     window = NULL,
                                     burnin=0, TOT=1000, k.max=6,
                                     stat=NULL, mu.alfa=NULL,
                                     mu.beta=NULL,
                                     s1=NULL, s2=NULL,  init.mu,
                                     init.sigma.2,
                                     init.beta, prob.k=NULL,
                                     sigma.tau.mu, sigma.tau.sigma.2, sigma.tau.beta,
                                     tau.split.mu, tau.split.beta,
                                     start.k, RJ=RJ, 
                                     NC, deltaT,
                                     temp.dir) {
    ##Hyperparameters
  if(is.null(mu.alfa)) mu.alfa <- median(y)
  if(is.null(mu.beta)) mu.beta <- diff(range(y))
  if(is.null(s1)) s1 <- mu.beta
  if(is.null(s2)) s2 <- mu.beta
  if(is.null(maxVar)) maxVar <- diff(range(y))^2
  
  ## Prior over the number of states  ####################

  if (is.null(prob.k)) {
    prob.k <- rep(1/k.max, k.max)
  }
  pb <- c(1,rep(1/2, k.max-2),0)
  ps <- c(1, rep(1/2, k.max-2),0)

  if (is.null(stat)) {
    for(i in 1:k.max)  stat <- c(stat, rep(1/i, i))
  }
  if(is.null(sigma.tau.mu) | is.null(sigma.tau.sigma.2) | is.null(sigma.tau.beta)
   | is.null(tau.split.mu)) {
      cat("Searching jump parameters...\n")
      if(length(y) != (length(x) + 1))
          stop("Length of vector of distances and data vector are different in get jump!")
      
      params <- get.jump(y=y, x=x, k.max=k.max, Chrom=Chrom, model=model,
                         var.equal=var.equal, mV=maxVar, 
                         normal.reference=normal.reference,
                         window = window,
                         prob.k=prob.k, pb=pb, ps=ps,  
                         mu.alfa=mu.alfa, mu.beta=mu.beta, stat=stat,
                         temp.dir = temp.dir)
      if (is.null(sigma.tau.mu)) sigma.tau.mu <- params$sigma.tau.mu
      if (is.null(sigma.tau.sigma.2)) sigma.tau.sigma.2 <- params$sigma.tau.sigma.2
      if (is.null(sigma.tau.beta)) sigma.tau.beta <- params$sigma.tau.beta
      if(is.null(tau.split.mu)) tau.split.mu <- mean(params$sigma.tau.mu) ^2
  }

  ## These checks have to live here: cannot be anywhere before
  ## prior params are set.
  if(k.max != length(sigma.tau.mu)) stop("k.max != length(sigma.tau.mu)")
  if(k.max != length(sigma.tau.sigma.2)) stop("k.max != length(sigma.tau.sigma.2)")
  if(k.max != length(sigma.tau.beta)) stop("k.max != length(sigma.tau.beta)")
  
  cat("    Starting Reversible Jump\n")
  if(length(y) != (length(x) + 1))
      stop("Length of vector of distances and data vector are different in metropolis!")
  res <- MetropolisSweep.C(y=y, x=x, k.max=k.max, Chrom=Chrom,
                           model=model, var.equal=var.equal,
                           maxVar=maxVar, 
                           normal.reference=normal.reference,
                           window = window,
                           burnin=burnin, TOT=TOT,
                           prob.k=prob.k, pb=pb, ps=ps,
                           mu.alfa=mu.alfa, mu.beta=mu.beta,
                           s1=s1, s2=s2,  init.mu=init.mu,
                           init.sigma.2=init.sigma.2,
                           init.beta=init.beta, 
                           sigma.tau.mu=sigma.tau.mu,
                           sigma.tau.sigma.2=sigma.tau.sigma.2,
                           sigma.tau.beta=sigma.tau.beta,
                           tau.split.mu=tau.split.mu,
                           tau.split.beta=tau.split.beta, stat=stat,
                           start.k=start.k, RJ=RJ,
                           NC=NC, deltaT=deltaT,
                           write_seq = 1,
                           temp.dir = temp.dir)
  res
}



scale.bound.Dist <- function(x, max.dist, nrowy) {
    ## Distance vector: scale, bound, and miscell checks
    ## would allow to use a max.dist that differs by chromosome
    x.old <- x
    ## if all distances are the same, they should be zero
    if (is.null(x) |
        isTRUE(all.equal(min(x.old, na.rm=TRUE),
                         max(x.old, na.rm=TRUE)))) {
        x <- rep(0, nrowy) ## zz:??? a "rep"?
    }
    ## Scale x'2 to avoid overflow
    if (!is.null(max.dist)) {
        x <- x/max.dist
    }
    else {
        xmax <- max(x, na.rm=TRUE)
        if(xmax < 0)
            stop("Negative distances. There is probably a problem with the Pos you used")
        if(xmax != 0)
            x <- x/max(x, na.rm=TRUE)
    }
    ## bound x, and delete NAs
    x[x > 1] <- 1
    x[x < 0] <- 0
    ## in Genome model, there are NAs, which are never used, but we do not
    ## want C to complaint
    x[is.na(x)] <- -1
    return(x)
}



RJaCGH <- function(y, Chrom=NULL, Start=NULL, End=NULL, Pos=NULL,
                   Dist=NULL, probe.names=NULL, maxVar=NULL,
                   model="Genome", var.equal=TRUE, max.dist=NULL,
                   normal.reference=0, 
                   window = NULL,
                   burnin=10000, TOT=10000, k.max=6,
                   stat=NULL, mu.alfa=NULL, mu.beta=NULL,
                   s1=NULL, s2=NULL, init.mu=NULL, init.sigma.2=NULL,
                   init.beta=NULL,
                   prob.k=NULL, jump.parameters=list(),
                   start.k=NULL, RJ=TRUE, 
                   NC=1, deltaT=2,
                   delete_gzipped = .__DELETE_GZIPPED,
                   singleState = FALSE) {
    cl <- match.call()
    sigma.tau.mu <- jump.parameters$sigma.tau.mu
    sigma.tau.sigma.2 <- jump.parameters$sigma.tau.sigma.2
    sigma.tau.beta <- jump.parameters$sigma.tau.beta
    tau.split.mu <- jump.parameters$tau.split.mu
    tau.split.beta <- NULL

    ## one or more arrays?
    if (is.null(dim(y))) {
        narrays <- 1
        y <- matrix(y, ncol = 1)
        colnames(y) <- "array1"
    } else {
        narrays <- ncol(y)
        if (is.null(colnames(y))) {
            colnames(y) <- paste("array", 1:ncol(y))
        }
    }

    if (!is.null(Pos)) {
      if (length(Pos) != nrow(y)) {
        stop ("Positions and log-ratios must have the same length\n")
      }
    }

    if (!is.null(Chrom)) {
      if (length(Chrom) != nrow(y)) {
        stop ("Chromosomes and log-ratios must have the same length\n")
      }
    }

    
    ## FIXME ZZ: verify no array name is the same as a name of
    ## any of the other components of the object. If it is, rename.
    ## set to null el colnames(y) and do el rep
    
    ############# Checks
    if (k.max < 2) stop("k.max must be 2 or more\n")
    if(!is.null(start.k)) {
        if(k.max < start.k)
            stop("start.k cannot be larger than k.max")
        if (length(init.mu) != start.k &
            length(init.sigma.2) != start.k &
            length(init.beta) != start.k * start.k) {
            stop ("lengths of starting values incorrect\n")
        }
    }
    check.no.float(Pos)
    check.no.float(Start)
    check.no.float(End)
    if (!is.null(start.k) & is.null(init.mu) & is.null(init.sigma.2) &
        is.null(init.beta)) {
        stop ("If start.k is selected, starting values must be provided\n")
    }
    
    ## Allow Chroms to have numbers or characters or a mix
    if(!is.null(Chrom)) {
        Chrom <- as.character(Chrom)
        uniqueChrom <- unique(Chrom)
    } else {
        Chrom <- rep(dummyChr, nrow(y))
        uniqueChrom <- dummyChr
    }
    
    ## Pos, Pos.rel, Distances, etc
##     Pos.rel <- NULL
    if(is.null(Pos) && (xor(is.null(Start), is.null(End)))) {
        stop("If Pos is NULL, BOTH Start and End must be provided")
    }
    if(is.null(Pos) && is.null(Start) && is.null(End)) {
        warning("Pos, Start, End are all NULL: ",
                "Pos created as vector of consecutive integers")
        Pos <- 1:nrow(y)
    }
    if(!is.null(Pos) && !is.null(Start) && !is.null(End)) {
        warning("Pos, Start, End are ALL non-NULL. ",
                "Only Pos will be used")
    }
    if(is.null(Pos) && !is.null(Start) && !is.null(End)) {
        Pos <- Start + (End - Start) / 2
        ## Sort data in case...
        if (!is.null(Dist)) {
            if(length(Dist) != (nrow(y) - 1))
               stop("Length of Dist should be number of observations minus 1")
            Dist <- c(Dist, NA)
        }
        ## In order to subset chromosome from Dist
        for (i in uniqueChrom) {
            if (any(Pos[Chrom==i] < 0)) {
                ids <- order(Pos[Chrom==i])
                if(!identical(ids, seq_along(ids)))
                    warning("Need to reorder data because of Pos, Start, End. ",
                            "In addition, a probably meaningless value of 0 might ",
                            "be assigned to one entry of Dist")
                y[Chrom==i,] <- y[Chrom==i,][ids,]
                Start[Chrom==i] <- Start[Chrom==i][ids]
                End[Chrom==i] <- End[Chrom==i][ids]
                Pos[Chrom==i] <- Pos[Chrom==i][ids]
                if (!is.null(probe.names)) {
                    probe.names[Chrom==i] <- probe.names[Chrom==i][ids]
                }
                if (!is.null(Dist)) {
                    Dist[Chrom==i] <- Dist[Chrom==i][ids]
                }
            }
        }
        if (!is.null(Dist)) {
            ## Take out last value and check last NA has not moved
            ## As we re-order Dist, the last position (NA)
            ## could have moved. We ensure Dist is the correct length,
            ## and set to 0 all NA.
            Dist <- Dist[-length(Dist)]
            Dist[is.na(Dist)] <- 0
        }
        if (is.null(Dist)) {
            ## Posible definition of distance between probes
            Dist <- Start[-1] - End[-length(End)]
        }
    }
    Pos.rel <- Pos
    if (!is.null(Chrom) && any(diff(Pos)<0)) {
        ## zz: FIXME: but Pos.rel only saved rarely!!!
        ## save relative positions for pREC_A, plots, etc.
##         Pos.rel <- Pos
        last.Pos <- 0
        for (i in 2:length(Pos)) {
            if (Chrom[i] != Chrom[i-1]) {
                last.Pos <- Pos[i-1]
            }
            Pos[i] <- Pos[i] + last.Pos
        }
    }
    chrom.Dist <- diff(Pos) - 1
    for (i in uniqueChrom) 
        chrom.Dist[Chrom==i][length(chrom.Dist[Chrom==i])] <- NA
    if(any(na.omit(chrom.Dist) < -1))
        stop("chrom.Dist < - 1; the Pos values used are not valid")

    if (!is.null(Dist)) chrom.Dist <- Dist
    x.for.model <- chrom.Dist

    if(model != "Chrom") {
        x.for.model <- scale.bound.Dist(x.for.model, max.dist, nrow(y))
    }
    if(model == "Chrom") {
        for(chr in seq_along(uniqueChrom)) {
            x.for.model[Chrom == chr] <-
                scale.bound.Dist(x.for.model[Chrom == chr], max.dist,
                                 sum(Chrom == chr))
        }
    }
    temp.dir <- tempdir2()
    out <- list()
    for(arr in colnames(y)) {
        cat(" Doing Array ", arr, "\n")
        ## This big loop we do not want to duplicate code for either
        ## iterating over Chroms (model == Chrom) or not. We always
        ## iterate, though might iterate only once if model not chrom. If
        ## model = "Chrom", we make as many fits as Chromosomes,
        ## otherwise, a single fit.

        ## If same model, this is fixed
        if(model != "Chrom") {
            yyC <- y[, arr]
            ChromC <- Chrom
            chrom.DistC <- x.for.model[-length(x.for.model)]
        }
        if(model == "Chrom") {
            iterChrom <- uniqueChrom
            tmp.viterbi <- list()
        } else iterChrom <- 1
        
        for(chr in iterChrom) { ## loop over chroms, if needed
            ## if Chrom, need to alter arguments
            if(model == "Chrom") {
                yyC <- y[Chrom == chr, arr]
                ChromC <- Chrom[Chrom == chr]
                chrom.DistC <- x.for.model[Chrom == chr]
                chrom.DistC <- chrom.DistC[-length(chrom.DistC)]
                cat("    Doing Chromosome", chr, "\n")
            }
            ## FIXME: code for filenames of viterbi could come here?
            res <- RJMCMC.NH.HMM.Metropolis(y=yyC, Chrom=ChromC,
                                            x=chrom.DistC,
                                            var.equal=var.equal,
                                            max.dist=max.dist, maxVar=maxVar, 
                                            normal.reference=normal.reference,
                                            window = window,
                                            burnin=burnin, TOT=TOT,
                                            k.max=k.max, stat=stat,
                                            mu.alfa=mu.alfa, mu.beta=mu.beta,
                                            s1=s1, s2=s2,  init.mu=init.mu,
                                            init.sigma.2=init.sigma.2,
                                            init.beta=init.beta, prob.k=prob.k,
                                            sigma.tau.mu=sigma.tau.mu,
                                            sigma.tau.sigma.2=sigma.tau.sigma.2,
                                            sigma.tau.beta=sigma.tau.beta,
                                            model=model,
                                            tau.split.mu=tau.split.mu,
                                            tau.split.beta=tau.split.beta,
                                            start.k=start.k, RJ=RJ,
                                            NC=NC, deltaT=deltaT,
                                            temp.dir = temp.dir)
            ret.list <- list(normal.reference = normal.reference,
                             window = window,
                             singleState = singleState,
                             prob.s = res$prob.s,
                             prob.b = res$prob.b,
                             prob.d = res$prob.d,
                             prob.c = res$prob.c,
                             prob.e = res$prob.e,
                             k = res$k,
                             y = res$y,
                             sdY = res$sdY,
                             lengthY = res$lengthY,
                             fit.k = res$fit.k
                             )
            if(model == "Chrom") {
                ## tmp.viterbi is used to force viterbi to be the last
                ## component in object
                tmp.viterbi[[chr]] <- get.viterbi.seq2(res)[[1]]
###                 out[[arr]]$viterbi[[chr]] <- get.viterbi.seq2(res)
                out[[arr]][[chr]] <- ret.list
                class(out[[arr]][[chr]]) <- "iRJaCGH"
            }
            else
                out[[arr]] <- ret.list
        } ## end loop over chroms
        
        if(model != "Chrom") {
            out[[arr]]$viterbi <- get.viterbi.seq2(res)
            class(out[[arr]]) <- c("RJaCGH.Genome", "iRJaCGH")
        }
        else {
            out[[arr]]$viterbi <- tmp.viterbi
            class(out[[arr]]) <- "RJaCGH.Chrom"
        }
    } ## end loop over arrays
    
    out$array.names <- colnames(y)
    out$Start <- Start
    out$End <- End
    out$Pos <- Pos
    out$Pos.rel <- Pos.rel
    out$Dist <- Dist
    out$Chrom <- Chrom
    out$Dist.for.model <- x.for.model
    out$probe.names <- probe.names
    out$call <- cl
    out$temp.dir <- temp.dir
##    out$Y <- y
    class(out) <- "RJaCGH"
    out
}

tempdir.RJaCGH <- function(obj) {
    cat("\n Directory for files associated to this object is ",
        obj$temp.dir, "\n")
}


check.no.float <- function(x) {
    ### Checks a vector is made up of integers
    our.float <- function(x) {
        if(is.integer(x))
            return(FALSE)
        else if(!is.null(x))
            !(isTRUE(all.equal(round(x), x)))
        else
            FALSE
    }
    if(our.float(x))
        stop(paste(deparse(substitute(x)),
                   "seems to cointain non-integer values"),
             call. = FALSE)
}


summary.iRJaCGH <- function(object, array=NULL, Chrom=NULL,
                             k=NULL, point.estimator="median",
                             quantiles=NULL, ...) {
  res <- list()
  if (is.null(k)) {
    k <- as.numeric(names(which.max(table(object$k))))
  }
  res$stat <- object$fit.k[[k]]$stat
  if (point.estimator=="mode") {
      dens <- apply(readRO2(object$fit.k[[k]]$mu), 2, density, bw="nrd0")
      res$mu <- unlist(lapply(dens, function(x) x$x[which.max(x$y)]))
      dim(res$mu) <- c(k, 1)
      colnames(res$mu) <- "mode"
      dens <- apply(readRO2(object$fit.k[[k]]$sigma.2), 2, density, bw="nrd0")
      res$sigma.2 <- unlist(lapply(dens, function(x)
                                   x$x[which.max(x$y)]))
      dim(res$sigma.2) <- c(k, 1)
      colnames(res$sigma.2) <- "mode"
      if(!(.__DELETE_BETA)) {
      dens <- apply(readRO2(object$fit.k[[k]]$beta), c(1,2), density, bw="nrd0")
      res$beta <- unlist(lapply(dens, function(x) x$x[which.max(x$y)]))
      res$beta <- matrix(res$beta, k)
      diag(res$beta) <- 0
      } else {
          res$beta <- NA
      }
    }
  else{

    if (is.null(quantiles)) {
      quantiles <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    }
    res$mu <- apply(matrix(readRO2(object$fit.k[[k]]$mu), ncol=k), 2, quantile,
                    probs=quantiles)
    res$mu <- t(res$mu)
    dim(res$mu) <- c(k, length(quantiles))
    colnames(res$mu) <- paste(round(100*quantiles), "%", sep="")
    res$sigma.2 <- apply(matrix(readRO2(object$fit.k[[k]]$sigma.2), ncol=k), 2,
                         quantile, probs=quantiles)
    res$sigma.2 <- t(res$sigma.2)
    dim(res$sigma.2) <- c(k, length(quantiles))
    colnames(res$sigma.2) <- paste(round(100*quantiles), "%", sep="")
    if(!(.__DELETE_BETA)) {
        res$beta <- apply(readRO2(object$fit.k[[k]]$beta), c(1,2), point.estimator)
    } else {
        res$beta <- NA
    }
  }
##  browser()

  if (!is.null(object$fit.k[[k]]$state.labels)) {
      rownames(res$mu) <- rownames(object$fit.k[[k]]$state.labels)
      rownames(res$sigma.2) <- rownames(object$fit.k[[k]]$state.labels)
      if(!(.__DELETE_BETA)) {
          rownames(res$beta) <- rownames(object$fit.k[[k]]$state.labels)
          colnames(res$beta) <- rownames(object$fit.k[[k]]$state.labels)
      }
      names(res$stat) <- rownames(object$fit.k[[k]]$state.labels)
  }
  res$k <- table(object$k)
  class(res) <- "summary.iRJaCGH"
  res
}


summary.RJaCGH.Chrom <- function(object, array=NULL, Chrom=NULL,
                                  k=NULL, point.estimator="median",
                                  quantiles=NULL, ...) {

  res <- list()
  if(length(Chrom) == 1) {
    res <- summary(object[[Chrom]], k=k, 
                   point.estimator=point.estimator, quantiles=quantiles)
  }
  else {
    for (chr in Chrom) {
      k <- as.numeric(names(which.max(table(object[[chr]]$k))))
      res[[chr]] <- summary.iRJaCGH(object[[chr]], array=array, Chrom=Chrom, k=k, 
                          point.estimator=point.estimator,
                          quantiles=quantiles)
    }
    names(res) <- Chrom
    class(res) <- "summary.RJaCGH.Chrom"
  }
  res
}

summary.RJaCGH.Genome <- function(object, array=NULL,
                                   Chrom=NULL, k=NULL,
                                   point.estimator="median",
                                   quantiles=NULL, ...) {
  res <- summary.iRJaCGH(object, array=array,
                          Chrom=Chrom, k=k, point.estimator, quantiles)
  res
}

summary.RJaCGH <- function(object, array=NULL, Chrom=NULL,
                            k=NULL, point.estimator="median",
                                 quantiles=NULL, ...) {
  res <- list()
  if (is.null(array)) {
    array <- object$array.names
  }
  if (is.null(Chrom)) {
    Chrom <- unique(object$Chrom)
  }
  for (i in array) {
    res[[i]] <- summary(object[[i]], array=array, Chrom=Chrom,
                        k=k, point.estimator=point.estimator,
                        quantiles=quantiles)
  }
  class(res) <- "summary.RJaCGH"
  res
}

print.summary.iRJaCGH <- function(x, ...) {
  cat("\nDistribution of the number of hidden states:\n")
  print(round(prop.table(x$k), 3), ...)
  cat("\nModel with ", length(x$stat), " states:\n")
  cat("\nDistribution of the posterior means of hidden states:\n")
  print(round(x$mu, 3), ...)
  cat("\nDistribution of the posterior variances of hidden states:\n")
  print(round(x$sigma.2, 3), ...)
  cat("\nParameters of the transition functions:\n")
  print(round(x$beta, 3), ...)
}

print.summary.RJaCGH.Genome <- function(x, ...) {
  print.summary.iRJaCGH(x, ...)
}

print.summary.RJaCGH.Chrom <- function(x, ...) {
  for(chr in names(x)) {
      cat("\nSummary for chromosome ", chr, ":\n")
      print(x[[chr]], ...)
      cat("\n---------------------------------------------\n")
  }
}

print.summary.RJaCGH <- function(x, ...) {
  for(arr in 1:length(x)) {
    cat("\nSummary for ARRAY ", names(x)[arr], ":\n")
    print(x[[arr]], ...)
    cat("\n================================================\n")
  }
}

smoothMeans <- function(obj, array=NULL, Chrom=NULL, k=NULL) {
  UseMethod("smoothMeans")
}

smoothMeans.iRJaCGH <- function(obj, array=NULL, Chrom=NULL, k=NULL) {
  n <- obj$lengthY # n <-  length(readRO2(obj$y))
  sumfit <- summary(obj$k)
  if (!is.null(k)) {
    if (sumfit[k] == 0) stop("No observations in that model\n")
    probs <- 1
  }
  else {
    k <- as.numeric(names(sumfit))
    probs <- prop.table(sumfit)
  }
  res <- matrix(0, n, length(k))
  for (i in 1:length(k)) {
    if(sumfit[i] > 0) {
##      res[,i] <- apply(obj$fit.k[[i]]$prob.states *
      prob.states <- readRO2(obj$fit.k[[k[i]]]$prob.states)
      mu <- readRO2(obj$fit.k[[k[i]]]$mu)
      res[,i] <- apply(prob.states *
                       matrix(rep(apply(mu, 2, mean),
                                  n), n, byrow=TRUE), 1, sum)
    }
    else {
      res[,i] <- 0
    }
  }
  res <- res %*% probs
  res
}

smoothMeans.RJaCGH.Genome <- function(obj, array=NULL, Chrom=NULL, k=NULL) {
  res <- smoothMeans.iRJaCGH(obj, array=array, Chrom=Chrom, k=NULL)
  res
}

smoothMeans.RJaCGH.Chrom <- function(obj, array=NULL,
                                      Chrom=NULL, k=NULL) {
  res <- list()
  for (chr in Chrom) {
    res[[chr]] <- smoothMeans.iRJaCGH(obj[[chr]], array=array,
                                       Chrom=Chrom, k=k)
  }
  res <- do.call("c", res)
  res
}

smoothMeans.RJaCGH <- function(obj, array=NULL, Chrom=NULL, k=NULL) {
  if (is.null(array)) {
    array <- obj$array.names
  }
  if (is.null(Chrom)) {
    Chrom <- unique(obj$Chrom)
  }
  if (inherits(obj[[array[1]]], 'RJaCGH.Genome')) {
    N <- obj[[array[1]]]$lengthY
  }
  else if (inherits(obj[[array[1]]], 'RJaCGH.Chrom')) {
    N <- 0
    for (i in Chrom) {
      N <- N + obj[[array[1]]][[i]]$lengthY
    }
  }
  
  res <- matrix(NA, nrow = N,
                ncol = length(array))
  for (i in seq_along(array)) {
    res[, i] <- smoothMeans(obj[[i]], Chrom=Chrom, k = k)
  }
  res
}

states <- function(obj, array=NULL, Chrom=NULL, k=NULL) {
  UseMethod("states")
}
## Now does the same with joint probs.
states.iRJaCGH <- function(obj, array=NULL, Chrom=NULL, k=NULL) {
  res <- NULL
  if (is.null(k)) {
    k <- as.numeric(names(which.max(table(obj$k))))
  }
  ##if (nrow(obj[[k]]$mu)==0) stop ("No observations in that HMM\n")
##  res$states <- apply(obj$fit.k[[k]]$prob.states, 1, which.max)
  
  res$states <- apply(readRO2(obj$fit.k[[k]]$prob.states), 1, which.max) ## This could easily be done in C?
  res$states <- factor(res$states, levels=1:k)
  res$prob.states <- readRO2(obj$fit.k[[k]]$prob.states)

    ## Region of normal, gain and loss
  if (is.null(obj$fit.k[[k]]$state.labels)) {
    ref <- as.numeric(names(which.max(table(res$states))))
    colnames(res$prob.states) <- 1:k
    colnames(res$prob.states)[ref] <- "Normal"
    levels(res$states)[ref] <- "Normal"
    if (ref < k) {
      colnames(res$prob.states)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
      levels(res$states)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
    }
    if (ref > 1) {
      colnames(res$prob.states)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
      levels(res$states)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
    }
  }
  else {
    idnames <- apply(obj$fit.k[[k]]$state.labels, 1, which.max)
    all.names <- colnames(obj$fit.k[[k]]$state.labels)[idnames]
    colnames(res$prob.states) <- all.names
    levels(res$states) <- all.names
  }
  res <- list(states=res$states, prob.states=res$prob.states)
  res
}

states.RJaCGH.Chrom <- function(obj, array=NULL, Chrom =NULL, k=NULL) {
  res <- list()
  for (i in Chrom) {
    if (is.null(k)) {
      k <- as.numeric(names(which.max(table(obj[[i]]$k))))
    }
    res[[i]] <- states(obj[[i]], array=array, Chrom=Chrom, k=k)
  }
  res
}

states.RJaCGH.Genome <- function(obj, array=NULL, Chrom = NULL, k=NULL) {
 states.iRJaCGH(obj, array=array, Chrom=Chrom, k=k)
}

states.RJaCGH <- function(obj, array=NULL, Chrom = NULL, k=NULL) {
  res <- list()
  if (is.null(array)) {
    array <- obj$array.names
  }
  if (is.null(Chrom)) {
    Chrom <- unique(obj$Chrom) 
  }
  for (i in array) {
    res[[i]] <- states(obj[[i]], array=array, Chrom = Chrom, k=k)
  }
  res
}

modelAveraging <- function(obj, array=NULL, Chrom=NULL) {
  UseMethod("modelAveraging")
}

modelAveraging.iRJaCGH <-function(obj, array=NULL, Chrom = NULL) {
  res <- list()
  n <- obj$lengthY # n <- length(readRO2(obj$y))
  probs <- prop.table(table(obj$k))
  Loss <- rep(0, n)
  Normal <- rep(0, n)
  Gain <- rep(0, n)
  
  k.max <- max(as.numeric(as.character(obj$k)))
  for (i in 1:k.max) {

    if (probs[i] > 0) {
      objSummary <- states.iRJaCGH(obj, k=i)
      prob.hiddenStates <- obj$fit.k[[i]]$state.labels
      prob.states <- objSummary$prob.states
      prob.states <- prob.states %*% prob.hiddenStates
      Loss <- Loss + probs[i] * prob.states[,1]
      Normal <- Normal + probs[i] * prob.states[,2]
      Gain <- Gain + probs[i] * prob.states[,3]
    }
  }
  
  res$prob.states <- cbind(Loss, Normal, Gain)
  colnames(res$prob.states) <- c("Loss", "Normal", "Gain")
  res$states <- apply(res$prob.states, 1, function(x) names(which.max(x)))
  res$states <- ordered(res$states, levels=c("Loss", "Normal", "Gain"))
  if (!is.null(obj$probe.names)) {
    names(res$states) <- obj$probe.names
    rownames(res$prob.states) <- obj$probe.names
  }
  res <- list(states=res$states, prob.states=res$prob.states)


}

modelAveraging.RJaCGH.Chrom <-function(obj, array=NULL, Chrom = NULL) {
  res <- list()
  for (ch in Chrom) {
    res[[ch]] <- modelAveraging(obj[[ch]])
  }
  res
}

modelAveraging.RJaCGH.Genome <- function(obj, array=NULL, Chrom = NULL) {
  modelAveraging.iRJaCGH(obj=obj, array=array, Chrom=Chrom)
}

modelAveraging.RJaCGH <- function(obj, array=NULL, Chrom = NULL) {
  res <- list()
  if (is.null(array)) {
    array <- obj$array.names
  }
  if (is.null(Chrom)) {
    Chrom <- unique(obj$Chrom)
  }
  for (i in array) {
    res[[i]] <- modelAveraging(obj[[i]], array=array, Chrom = Chrom)
  }
  res
}






plot.RJaCGH.Genome <- function(x, array=NULL, k=NULL,  Dist.for.model=NULL, 
                               model.averaging=TRUE, cex=1,
                               smoother=FALSE, Pos=NULL,
                                Pos.rel = NULL, Chrom=NULL, ...)  {
  layout(matrix(c(1,3,5,2,4, 5),3,2), heights=c(0.25, 0.25, 0.5))
  opar <- par(cex.lab=cex*0.8,   cex.main=cex*0.9,
              cex.axis=cex*0.8,  mar=c(5,4,4,4) + 0.1)
  
  y <- readRO2(x$y)
    
    if (is.null(k)) {
        k <- as.numeric(names(which.max(table(x$k))))
    }
    
    if (model.averaging) {
        res <- modelAveraging(x)
        cols <- colnames(states(x, k)$prob.states)
    }
    else {
        res <- states(x, k=k)
        cols <- colnames(res$prob.states)
    }
    
    barplot(prop.table(table(x$k)), xlab="# of hidden states", ylab="probability",
            main="Prob. of number of hidden states", col="red", ylim=c(0,1))
    
    if(model.averaging) {
        names.colors <- colnames(states(x, k)$prob.states)
        col <- rep(1, k)
        col[grep("[G]", names.colors)] <- 2
        col[grep("[L]", names.colors)] <- 3
        
    } else {
        col <- rep(1, k)
        col[grep("[G]", cols)] <- 2
        col[grep("[L]", cols)] <- 3
    }
    mu <- readRO2(x$fit.k[[k]]$mu)
    sigma.2 <- readRO2(x$fit.k[[k]]$sigma.2)
    plot(density(mu[,1], bw=sd(mu[,1])), col=col[1], xlim=range(y),
         main="Posterior probability of mean of hidden states")
    if (k >1) for (i in 2:k) lines(density(mu[,i],
                                           bw=sd(mu[,i])), col=col[i])
    
    plot(density(sigma.2[,1],
                 bw=sd(sigma.2[,1]), from=0), col=col[1], main="Posterior probability of variance of hidden states")
    if (k >1) for (i in 2:k) lines(density(sigma.2[,i],
                                           bw=sd(sigma.2[,i])), col=col[i])
    
    summary.obj <- summary(x, k)
  plotQNH(q=-summary.obj$beta, beta=summary.obj$beta, x=Dist.for.model[!is.na(Dist.for.model)],
              main="Probability of permanence in the same hidden state", xlab="Distance", ylab="Prob.", col=col)
  internal.plot.obs.pred(y = y,
                           Pos = Pos,
                           Pos.rel = Pos.rel,
                           states = res$states,
                           Chrom = Chrom,
                           model.averaging = model.averaging,
                           prob.states = res$prob.states,
                           smoother = smoother,
                           cex = cex, k=k)
    
}

internal.plot.obs.pred <- function(y, Pos, Pos.rel, states, Chrom,
                                   model.averaging, prob.states,
                                   smoother, cex, k, prob=TRUE) {
    if (model.averaging)
        main.text <- "Prediction of copy gain/loss. Bayesian Model Averaging"
    else
        main.text <- paste("Prediction of copy gain/loss. Number of hidden states: ", k)
    ylim <- range(y)
    margin <- diff(ylim) / 4
    ylim[1] <- ylim[1] - margin - 1
    ylim[2] <- ylim[2] + margin
    col <- rep(1, length(y))
    col[as.numeric(states) %in% grep("[G]", levels(states))] <- 2
    col[as.numeric(states) %in% grep("[L]", levels(states))] <- 3
    pch <- rep(16, length(y))
    pch[as.numeric(states) %in% grep("[2-9]", levels(states))] <- 17
    
    if(!is.null(Chrom)) {
        ## Start of every Chromosome
        start.Chrom <- c(Pos[which(!duplicated(Chrom))],
                         Pos[length(Chrom)] + 1)
        xx <- rep(start.Chrom, rep(2, length(start.Chrom)))
        ## Small hack for the case of even number of chroms
        if (length(start.Chrom) %% 2 !=0) {
            xx <- c(xx, xx[2])
        }
        yy <- rep(c(ylim, rev(ylim)), length.out=length(xx))
        ##         Pos <- x$Pos
    } else {
        ##         Pos <- x$Pos
        if (!is.null(Pos.rel)) Pos <- Pos.rel
    }
    
    plot(y ~ Pos, pch=pch, col=col, ylim=ylim, ylab="Log2 ratio", xlab="Pos.Base",
         main=main.text, type="n")
    
    if(!is.null(Chrom)) {
        ## here, or points will not be visible
        polygon(xx, yy, col="gray85", xpd=TRUE)
        text(x = start.Chrom[-length(start.Chrom)] + diff(start.Chrom)/2,
             y = rep(ylim[2] , length(start.Chrom)), labels=unique(Chrom),
             cex=cex*0.6, pos=1)
    }
    points(y ~ Pos, pch=pch, col=col)
    abline(h=0, lty=2)
    if (smoother) {
        smoothed.means <- smoothMeans(y)
        lines(smoothed.means~Pos, col="orange")
    }
    if (prob) {
      par(new=TRUE)
      plot(apply(prob.states, 1, max) ~ Pos, yaxt="n", xlab="", ylab="",
           col="blue", ylim=c(0,5), type="s")
      abline(h=c(0, 0.25, 0.5, 0.75,1), lty=2, yaxt="n")
      axis(4, at=c(0, 0.5, 0.75, 1), labels=c(0, 0.5, 0.75, 1))
      mtext("Probability", side = 4, line = 2, cex = cex*0.6, at = 0.5)
    }
}


plot.RJaCGH.Chrom <- function(x, Chrom=NULL, array=NULL,
                              model.averaging=TRUE, cex=1, k=NULL,
                              smoother=FALSE, ...)  {
  if (is.null(array)) stop ("Must specify array\n")
  if (is.null(Chrom)) {
        states <- NULL
        prob.states <- NULL
        ## Always model.averaging if Chrom==NULL
        model.averaging <- TRUE
        res <- modelAveraging(x)
        ## FIXME: cambiara algo de esto cuando guardemos summaires
        ## en el objeto? Creo que no. Pero estar al loro!!!
        max.levels <- lapply(res[[array]], function(x) levels(x$states))
        max.levels <- names(table(unlist(max.levels)))
        states <- lapply(res[[array]], function(x) as.character(x$states))
        states <- factor(unlist(states), levels=max.levels)
        prob.states <- lapply(res[[array]], function(x) x$prob.states)
        prob.states <- do.call("rbind", prob.states)
        par(cex.lab=cex*0.8, cex.main=cex*0.9, cex.axis=cex*0.8)
        par(mar=c(5,4,4,4) + 0.1)

        ##         y <- NULL
        ##         for (chr in unique(x$Chrom)) {
        ##             y <- c(y, x[[chr]]$y)
        ##         }
        ncrom <- length(unique(x$Chrom))

        y <- numeric(length(x$Pos))
        id <- 1
        for (i in 1:ncrom) {
          tmp <- readRO2(x[[array]][[i]]$y)
          y[id:(id + length(tmp)-1)] <- tmp
          id <- id + length(tmp)
        }
        Pos <- x$Pos
        internal.plot.obs.pred(y = y, Pos = Pos, Pos.rel = NULL,
                               states = states, Chrom = x$Chrom,
                               model.averaging = model.averaging,
                               prob.states = prob.states,
                               smoother = smoother,
                               cex = cex, k=k)
  }
  else {
    plot.RJaCGH.Genome(x[[array]][[Chrom]], array=array,
                        Dist.for.model = x$Dist.for.model, 
                        model.averaging=model.averaging,
                        cex=cex, k=k, smoother = smoother,
                        Pos=x$Pos.rel[x$Chrom==Chrom],
                        Chrom=x$Chrom[x$Chrom==Chrom])
  }
}




## There could be a weights argument for array precision
plot.RJaCGH <- function(x, array=NULL, k=NULL, Chrom=NULL, show="average", weights=NULL,
                         model.averaging=TRUE, 
                         cex=1, smoother=FALSE,...)  {
  
  if (!is.null(array)) {
    if(inherits(x[[array]], 'RJaCGH.Genome'))
      plot.RJaCGH.Genome(x=x[[array]], k=k, Dist.for.model = x$Dist.for.model, 
                          model.averaging=model.averaging, cex=cex,
                          smoother=smoother, Pos=x$Pos, Chrom=x$Chrom,...)
    if(inherits(x[[array]], 'RJaCGH.Chrom'))
      plot.RJaCGH.Chrom(x=x, array=array, k=k, Chrom=Chrom, model.averaging=
                          model.averaging, cex=cex, smoother=smoother,...)

  }
  else {
  

    if (show!="average" & show!="frequency") {
        stop("'show' must be either 'average' or 'frequency'\n")
    }
    uChrom <- unique(x$Chrom)
    array.names <- x$array.names
    par(mfrow=c(1,1))
    
    ##To make apply only on arrays we take out element 'array.names'

    ## FIXME!! Ojo, no, mejor hace apply solo en los elementos apropiados
    ## (1:numero arrays) y asi evitamos modificar x.
    
    # x$array.names <- NULL
    if (is.null(weights)) weights <- rep(1/length(array.names), length(array.names))
    weights <- weights / sum(weights)
    model.averaging <- TRUE
    res <- modelAveraging(x)
    y <- matrix(0, length(x$Pos), length(x$array.names))
    prob.states <- matrix(0, length(x$Pos), 3)

    if(inherits(x[[array.names[1]]], 'RJaCGH.Genome')) {
      res <- lapply(res, function(x) x$prob.states)
      for (i in 1:length(array.names)) {
        prob.states <- prob.states + res[[i]] * weights[i]
        y[,i] <- readRO2(x[[array.names[i]]]$y)
      }

      states <- apply(prob.states, 1, which.max)
      states <- factor(states, levels=1:3, labels=
                       c("Loss", "Normal", "Gain"))
    }
    else if(inherits(x[[array.names[[1]]]], 'RJaCGH.Chrom')) {
      ncrom <- length(unique(x$Chrom))
      res <- lapply(res, function(x) lapply(x, function(y) y$prob.states))
      for (i in 1:length(array.names)) {
        id <- 1
        prob.states <- prob.states + do.call("rbind",res[[i]]) *
          weights[i]
        for (j in 1:ncrom) {
          tmp <- readRO2(x[[array.names[i]]][[j]]$y)
          y[id:(id + length(tmp)-1),i] <- tmp
          id <- id + length(tmp)
        }
        id <- 1
      }
      states <- apply(prob.states, 1, which.max)
      states <- factor(states, levels=1:3, labels=
                       c("Loss", "Normal", "Gain"))
    }
    
    if (show=="average") {
      y <- y %*% weights
      internal.plot.obs.pred(y = y, Pos = x$Pos, Pos.rel = NULL,
                             states = states, Chrom = x$Chrom,
                             model.averaging = model.averaging,
                             prob.states = prob.states,
                             smoother = smoother,
                             cex = cex, k=k)
    }
    else {
      code <- apply(prob.states[,-2], 1, which.max)
      max.prob.states <- rep(0, nrow(prob.states))
      max.prob.states[code == 1] <- -prob.states[code == 1,1]
      max.prob.states[code == 2] <- prob.states[code == 2,3]
      internal.plot.obs.pred(y = max.prob.states, Pos = x$Pos,
                             Pos.rel = NULL,
                             states = states, Chrom = x$Chrom,
                             model.averaging = model.averaging,
                             prob.states = prob.states,
                             smoother = smoother,
                             cex = cex, k=k, prob=FALSE)
    }
  }
}

trace.plot <- function(x, k=NULL, array=NULL, Chrom=NULL, main.text=NULL) {

  if (is.null(array)) {
    stop ("'array' must be specified\n")
  }
  if(inherits(x[[array]], 'RJaCGH.Genome')) {
    mod <- x[[array]]$k
    if (is.null(k)) {
      k <- as.numeric(names(which.max(table(x[[array]]$k))))
    }
    mu <- readRO2(x[[array]]$fit.k[[k]]$mu)
    sigma.2 <- readRO2(x[[array]]$fit.k[[k]]$sigma.2)
    beta <- readRO2(x[[array]]$fit.k[[k]]$beta)
  }
  if(inherits(x[[array]], 'RJaCGH.Chrom')) {
    if (is.null(Chrom)) {
      stop ("'Chrom' must be specified\n")
    }
    mod <- x[[array]][[Chrom]]$k
    if (is.null(k)) {
      k <- as.numeric(names(which.max(table(x[[array]][[Chrom]]$k))))
    }
    mu <- readRO2(x[[array]][[Chrom]]$fit.k[[k]]$mu)
    sigma.2 <- readRO2(x[[array]][[Chrom]]$fit.k[[k]]$sigma.2)
    beta <- readRO2(x[[array]][[Chrom]]$fit.k[[k]]$beta)
  }
  ## Something must be made with colors
  par(mfrow=c(2, 2))
  matplot(as.numeric(as.character(mod)), pch=16, cex=0.2, main="Trace plot of number of states",
          xlab="iteration", ylab="Number of states", type="l")
  matplot(mu, type="l", main="Trace plot of means", xlab="iteration",
          ylab="Mean of states")
  matplot(sigma.2, type="l", main="Trace plot of variance", xlab="iteration",
          ylab="Variance of the states")
  if (k >1) {
    if(!(.__DELETE_BETA)) {
      matplot(t(beta[1,,]), type="n", main="Trace plot of beta", xlab="iteration",
              ylab="Beta")
      for (i in 1:k)
        matplot(t(beta[i,,]), type="l", add=TRUE)
    }
  }
  if (is.null(main.text)) {
    main.text <- paste("Whole genome.", k, "hidden states")
  }
  else {
    main.text <- paste(main.text, k, "hidden states")
  }
  mtext(main.text, outer=TRUE)
}


##############################
## Adapt parameters intra model
get.jump <- function(y, x, Chrom, model, k.max=6, normal.reference,
                     window,
                     ## normal.ref.percentile,
                     var.equal=TRUE,
                     mV=diff(range(y)), 
                     prob.k=NULL, pb=NULL, ps=NULL,
                     s1=NULL, s2=NULL,
                     init.mu=NULL,
                     init.sigma.2=NULL,
                     init.beta=NULL,
                     mu.alfa=NULL, mu.beta=NULL, stat,
                     increment=2, max.tries=10,
                     temp.dir) {
 
  TOL <- 0.01
  increment.mu <- increment
  increment.sigma.2 <- increment
  increment.beta <- increment
  optim.mu <- rep(0, k.max)
  optim.sigma.2 <- rep(0, k.max)
  optim.beta <- rep(0, k.max)
  for (k in 2:k.max) {
###       set.seed(SUPERSEED + k)
###       cat("\n            RUNIF in get.jump ", runif(1), "\n")

    ##   cat("max", sigma.tau.mu.max, "\n")
    ##   cat("mu", sigma.tau.mu, "\n")
    ##   cat("min", sigma.tau.mu.min, "\n")
    p.mu <- 0
    p.sigma.2 <- 0
    p.beta <- 0
    
    prob.lim <- c(0.2, 0.5)

    sigma.tau.mu <- 0.005
    sigma.tau.sigma.2 <- 0.005
    sigma.tau.beta <- 0.01

    tries <- 1

    ## Overdispersed init
    init.mu <- runif(k, -2, 2)
    if (!var.equal) {
        init.sigma.2 <- runif(k, 0, mV)
    }
    else {
        init.sigma.2 <- rep(runif(1, 0, mV), k)
    }
    init.beta <- matrix(rgamma(k * k, 1, 1), k)
    diag(init.beta) <- 0
    
    
    ## Check if min == max
  while ((!p.mu | !p.sigma.2 | !p.beta) & tries < max.tries) {
###        set.seed(SUPERSEED + k + 8)
###       cat("\n            RUNIF in get.jump before metropolis", runif(1), "\n")

    fit <- MetropolisSweep.C(y=y, x=x, k.max=k.max, Chrom=Chrom,
                             model=model, var.equal=var.equal,
                             maxVar=mV, 
                             burnin=0, TOT=1000,
                             normal.reference=normal.reference,
                             window = window,
##                             normal.ref.percentile=normal.ref.percentile,
                             prob.k=prob.k, pb=pb, ps=ps,
                             mu.alfa=mu.alfa, mu.beta=mu.beta,
                             s1=s1, s2=s2,
                             init.mu=init.mu,
                             init.sigma.2=init.sigma.2,
                             init.beta=init.beta,
                             sigma.tau.mu=rep(sigma.tau.mu,k.max),
                             sigma.tau.sigma.2=rep(sigma.tau.sigma.2, k.max),
                             sigma.tau.beta=rep(sigma.tau.beta, k.max),
                             stat=stat, RJ=FALSE, start.k=k, NC=1, deltaT=2,
                             write_seq = 0,
                             delete_gzipped = TRUE,
                             temp.dir = temp.dir)
    prob.mu <- fit$fit.k[[k]]$prob.mu
    prob.sigma.2 <- fit$fit.k[[k]]$prob.sigma.2
    prob.beta <- fit$fit.k[[k]]$prob.beta
##     cat("par=", sigma.tau.mu, " prob=", prob.mu, "\n")
##     cat("le voy a multiplicar por ", increment.mu, "\n")

    ## Check mu prob. of acceptance
    if (prob.mu < prob.lim[1] & !p.mu) {
      increment.mu <- increment * (prob.lim[1] / prob.mu)
      sigma.tau.mu <- sigma.tau.mu / increment.mu
    }
    if (prob.mu > prob.lim[2] | !p.mu) {
      increment.mu <- increment * (prob.mu / prob.lim[2])
      sigma.tau.mu <- sigma.tau.mu * increment.mu
    }
    
    if (prob.mu > prob.lim[1] & prob.mu < prob.lim[2] & !p.mu) {
      p.mu <- 1
    }
    ## Check sigma.2 prob. of acceptance
    if (prob.sigma.2 < prob.lim[1] & !p.sigma.2) {
      increment.sigma.2 <- increment * (prob.lim[1] / prob.sigma.2)
      sigma.tau.sigma.2 <- sigma.tau.sigma.2  / increment.sigma.2
    }
    if (prob.sigma.2 > prob.lim[2] | !p.sigma.2) {
      increment.sigma.2 <- increment * (prob.sigma.2 / prob.lim[2])
      sigma.tau.sigma.2 <- sigma.tau.sigma.2 * increment.sigma.2
    }
    
  if (prob.sigma.2 > prob.lim[1] & prob.sigma.2 < prob.lim[2] & !p.sigma.2) {
    p.sigma.2 <- 1
  }
  ## Check beta prob. of acceptance
  if (prob.beta < prob.lim[1] & !p.beta) {
    increment.beta <- increment * (prob.lim[1] / prob.beta)
    sigma.tau.beta <- sigma.tau.beta / increment
  }
  if (prob.beta > prob.lim[2] | !p.beta) {
    increment.beta <- increment * (prob.beta / prob.lim[2])
    sigma.tau.beta <- sigma.tau.beta * increment.beta
  }

  if (prob.beta > prob.lim[1] & prob.beta < prob.lim[2] & !p.beta) {
    p.beta <- 1
  }
##   cat("Iteration=", i, "\n")
##   cat("max", sigma.tau.mu.max, "\n")
##   cat("mu", sigma.tau.mu, "\n")
##   cat("min", sigma.tau.mu.min, "\n")
    tries <- tries + 1

  }
  optim.mu[k] <- sigma.tau.mu
  optim.sigma.2[k] <- sigma.tau.sigma.2
  optim.beta[k] <- sigma.tau.beta
}
  optim.mu[1] <- optim.mu[2]
  optim.sigma.2[1] <- optim.sigma.2[2]
  optim.beta[1] <- optim.beta[2]
  res <- list(sigma.tau.mu=abs(optim.mu), sigma.tau.sigma.2=abs(optim.sigma.2),
              sigma.tau.beta=abs(optim.beta))
  print(res)
  res
}

akaike <- function(logliks, param=NULL) {
  if (is.null(param)) {
    param <- 1:length(logliks)
    param <- 2 * param + param * (param-1)
  }
  AIC <- -2*logliks + 2*param
  AIC <- AIC - min(AIC)
  akaike <- exp(-0.5*AIC)
  akaike <- akaike / sum(akaike)
  akaike
}

genomePlot <- function(obj, array=NULL, weights=NULL, 
                       col=NULL, breakpoints=NULL, legend.pos=NULL,...) {

  if (!is.null(col) & !is.null(breakpoints)) {
    if(length(col) != length(breakpoints) + 1)
      stop("length(col) must be length(breakpoints + 1\n")
    if(length(breakpoints) < 2)
      stop("length(breakpoints) must be at least 2\n")
   }


  if (is.null(array)) {
    array <- obj$array.names
  }
  if (is.null(weights))
    weights <- rep(1/length(array), length(array))
  weights <- weights / sum(weights)
  res <- modelAveraging(obj)
  y <- matrix(0, length(obj$Pos), length(array))
  probs <- matrix(0, length(obj$Pos), 3)
  Pos <- obj$Pos.rel
  if(inherits(obj[[array[1]]], 'RJaCGH.Genome')) {
    res <- lapply(res, function(x) x$prob.states)
    for (i in 1:length(array)) {
      probs <- probs + res[[i]] * weights[i]
      y[,i] <- readRO2(obj[[array[i]]]$y)
      }
    
  }
  else if(inherits(obj[[array[[1]]]], 'RJaCGH.Chrom')) {
    ncrom <- length(unique(obj$Chrom))
    res <- lapply(res, function(x) lapply(x, function(y) y$prob.states))
      for (i in 1:length(array)) {
        id <- 1
        probs <- probs + do.call("rbind",res[[i]]) *
          weights[i]
        for (j in 1:ncrom) {
          tmp <- readRO2(obj[[array[i]]][[j]]$y)
          y[id:(id + length(tmp)-1),i] <- tmp
          id <- id + length(tmp)
        }
        id <- 1
      }
  }

  probs <- round(probs, 2)
  ## what if c(0, 0.5, 0.5) ? which.max=2
  states <- apply(probs, 1, which.max)
  ## Assign he normal probes to gain if y>0 or loss if y<0
  index <- 1*(y<0) + 3*(y>=0)
  states[states==2] <- index[states==2]
  index <- states
  colo <- rep(0, length(index))
  for(i in 1:length(index)) {
    colo[i] <- probs[i,index[i]]
  }
  colo[index==1] <- -colo[index==1]
  colo.round <- floor(colo*10) / 10
  colo.recoded <- rep(0, length(colo.round))
  
    if (is.null(col)) {
    col <- colors()
    col <- col[c(86, 50, 51, 24, 555, 552, 404)]
  }

  if(is.null(breakpoints)) {
    breakpoints <- c(-0.9, -0.7, -0.5, 0.5, 0.7, 0.9)
  }
  MidPoint <- floor(length(col)/2)
  colo.recoded[colo.round <= -breakpoints[1]] <- col[1]
  label.legend <- paste("P.Loss >= ", -breakpoints[1], sep="")
  if (length(breakpoints) > 2) {
    for (i in 2:MidPoint) {
      colo.recoded[colo.round > breakpoints[i-1] & colo.round <=
    breakpoints[i]] <- col[i]
      label.legend <- c(label.legend, paste(-breakpoints[i],
                                            " <= P.Loss < ", -breakpoints[i-1], sep=""))
    }
  }
  colo.recoded[colo.round > breakpoints[MidPoint] & colo.round <
    breakpoints[MidPoint +1]] <- col[MidPoint +1]
  label.legend <- c(label.legend, paste("P.Loss < ", -breakpoints[MidPoint],
                                        " or P.Gain < ", breakpoints[MidPoint+1], sep=""))
  if (length(breakpoints) > 2) {
    for (i in (MidPoint+2):length(breakpoints)) {
      colo.recoded[colo.round >=breakpoints[i-1] & colo.round
    < breakpoints[i]] <- col[i]
      label.legend <- c(label.legend, paste(breakpoints[i-1],
                                            " <= P.Gain < ", breakpoints[i], sep=""))
    }
  }

  colo.recoded[colo.round >=breakpoints[length(breakpoints)]] <-
  col[length(col)]
  label.legend <- c(label.legend, paste("P.Gain >= ",
                                        breakpoints[length(breakpoints)], sep=""))
  Chrom <- obj$Chrom
  n.chrom <- length(unique(Chrom))
  par(mar=c(0,0,0,2), oma=c(0, 4, 4, 4))
  layout(rbind(c(1, 2), c(3, 3)), widths=c(1,1), heights=c(9,2.5))
  xmax <- max(Pos)
  plot(0,0, type="n", xlim=c(0, xmax), ylim=c(1, ceiling(n.chrom/2)),
       axes=FALSE, ylab="Chromosome", xlab="",...)
  axis(side=2, at=c(1:ceiling(n.chrom/2)), labels=unique(Chrom)[1:ceiling(n.chrom/2)])

  chrom.count <- 1

  for(i in unique(Chrom)[1:ceiling(n.chrom/2)]) {
    lines(c(0, max(Pos[Chrom==i])), c(chrom.count,chrom.count))
    points(Pos[Chrom==i], chrom.count + y[Chrom==i]/ (2*max(abs(y))),
           pch=19,col=colo.recoded[Chrom==i], cex=0.5)
    chrom.count <- chrom.count + 1
  }
  plot(0,0, type="n", xlim=c(0, xmax), ylim=c(1, ceiling(n.chrom/2)),
       axes=FALSE, ylab="Chromosome", xlab="")
  axis(side=2, at=1:(n.chrom - ceiling(n.chrom/2)), labels=unique(Chrom)[(ceiling(n.chrom/2)+1):n.chrom])

  chrom.count <- 1
  for(i in unique(Chrom)[(ceiling(n.chrom/2) +1):n.chrom]) {

    lines(c(0, max(Pos[Chrom==i])), c(chrom.count,chrom.count))
    points(Pos[Chrom==i], chrom.count  + y[Chrom==i]/ (2*max(abs(y))),
           pch=16,col=colo.recoded[Chrom==i], cex=0.5)
    chrom.count <- chrom.count + 1
  }

  if(is.null(legend.pos)) {
    plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
    legend("bottom", legend=label.legend,
           col=col, pch=19, cex=0.9)
  }
  else {
    plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
    legend(x=legend.pos[1], y=legend.pos[2], legend=label.legend,
           col=col, pch=19, cex=0.9)
  }
}





relabelStates <- function(obj, normal.reference=0, window=NULL,
                          singleState = FALSE, array=NULL, Chrom=NULL) {
    ## singleState = TRUE assigns each state to a single state,
    ## so each state has a p = 1 of being something and
    ## a p = 0 for everything else. Emulates the old
    ## relabelStates.
  UseMethod("relabelStates")
}

relabel.core <- function(obj, normal.reference=0, window=NULL,
                         singleState = FALSE)  {
    if (is.null(window)) {
        window <- obj$sdY
    } else {
        window <- window * obj$sdY
    }
    k.max <- max(as.numeric(levels(obj$k)))
    model.probs <- prop.table(table(obj$k))
  state.labels.list <- list()
  for (i in 1:k.max) {
      if (model.probs[i] > 0) {
      if(!is.null(dim(obj$fit.k[[i]]$state.labels))) {
        obj.sum <- summary.iRJaCGH(obj, k=i, point.estimator="median",
                                  quantiles = 0.5)
      }
      else {
          obj.sum <- list()
          obj.sum$mu <- apply(readRO2(obj$fit.k[[i]]$mu), 2, median)
          obj.sum$sigma.2 <- apply(readRO2(obj$fit.k[[i]]$sigma.2), 2, median)
      }
      limits <- normal.reference + c(-1, 1) * window
      probs <- matrix(0, nrow=i, ncol=3)
      probs[,1] <- pnorm(limits[1], obj.sum$mu,
                         sqrt(obj.sum$sigma.2),
                         lower.tail=TRUE)
      probs[,3] <- pnorm(limits[2], obj.sum$mu,
                         sqrt(obj.sum$sigma.2),
                         lower.tail=FALSE)
      probs[,2] <- 1 - probs[,1] - probs[,3]
      state.labels.list[[i]] <- probs
      modal.state <- apply(probs, 1, which.max)
      
      if(singleState) {
          msarray <- cbind(seq_along(modal.state), modal.state)
          state.labels.list[[i]] <- matrix(0, nrow = length(modal.state),
                                          ncol = 3)
          state.labels.list[[i]][msarray] <- 1
      }
      colnames(state.labels.list[[i]]) <- c("Loss", "Normal", "Gain")
      ## rownames of state.labels
      rownames(state.labels.list[[i]]) <- letters[1:nrow(probs)]
      basicLabel <- c("Loss", "Normal", "Gain")
      for(ii in c(1, 2, 3)) {
          this.type <- which(modal.state == ii)
          l.this.type <- length(this.type)
          if(l.this.type == 0) next
          this.label <- basicLabel[ii]
          if(l.this.type == 1)
              rownames(state.labels.list[[i]])[this.type] <- this.label
          else
              rownames(state.labels.list[[i]])[this.type] <-
                  paste(this.label, 1:l.this.type, sep = "-")
      }      
    } else {## did not visit these states
        state.labels.list[[i]] <- matrix(rep(-9, 3 * i), ncol = 3)
        colnames(state.labels.list[[i]]) <- c("Loss", "Normal", "Gain")
    }
  }
  state.labels.list
}
                           
relabelStates.iRJaCGH <- function(obj, normal.reference=0, window=NULL,
                             singleState = FALSE,
                                   array=NULL, Chrom=NULL)  {
    sllist <- relabel.core(obj = obj,
                             normal.reference = normal.reference,
                             window = window,
                             singleState = singleState)
    k.max <- max(as.numeric(levels(obj$k)))
    for(i in (1:k.max)) {
        obj$fit.k[[i]]$state.labels <- sllist[[i]]
    }
    return(obj)
}

relabelStates.RJaCGH.Chrom <- function(obj, normal.reference=0, window=NULL,
                                        singleState = FALSE,
                                        array=NULL, Chrom=NULL)  {
  for(chr in Chrom) {
    obj[[chr]] <- relabelStates.iRJaCGH(obj[[chr]], normal.reference=
                                         normal.reference, window=window,
                                         singleState = singleState)
  }
  obj
}

relabelStates.RJaCGH.Genome <- function(obj, normal.reference=0, window=NULL,
                                    singleState = FALSE,
                                         array=NULL, Chrom=NULL)  {
    obj <- relabelStates.iRJaCGH(obj, normal.reference=
                                normal.reference, window=window,
                                singleState = singleState)
    obj
}


relabelStates.RJaCGH <- function(obj, normal.reference=0, window=NULL,
                                   singleState = FALSE,
                                  array=NULL, Chrom=NULL)  {
  if (is.null(array)) {
    array <- obj$array.names
  }
  if (is.null(Chrom)) {
    Chrom <- unique(obj$Chrom)
  }
  for (i in array) {
    obj[[i]] <- relabelStates(obj[[i]], normal.reference=
                              normal.reference, window=window,
                              singleState = singleState,
                              array=array, Chrom=Chrom)
  }
  obj
}



pREC_A <- function(obj, p, alteration = "Gain",
                   array.weights = NULL,
                   verbose = FALSE) {
    pREC(method = "pREC_A", obj = obj,
         p = p, alteration = alteration,
         array.weights = array.weights,
         verbose = verbose)
}


pREC_S <- function(obj, p, freq.array, alteration = "Gain",
                   verbose = FALSE) {
    pREC(method = "pREC_S", obj = obj,
         p = p, alteration = alteration,
         freq.array = freq.array,
         verbose = verbose,
         array.weights = NULL)
}





pREC <- function(method, obj, p, freq.array = NULL,
                 alteration = "Gain", array.weights = NULL,
                 verbose = FALSE) {

    if (method == "pREC_A") {
        method_prec <- 0
        freq.array <- -9
    } else if (method == "pREC_S") {
        method_prec <- 1
    } else {
        stop ("method can only be one of pREC_A or pREC_S")
    }
    
    if((method == "pREC_S") & is.null(freq.array))
        stop("pREC_S requires a value for freq.array")
    if(alteration == "Gain") {
        alteration.int <- 1
    } else if (alteration == "Loss") {
        alteration.int <- -1
    } else {
        stop("alteration can only be one of 'Gain' or 'Loss'")
    }

    array.names <- obj$array.names
    narrays <- length(obj$array.names)
    idchr <- unique(obj$Chrom)
    nchrom <- length(idchr)

    if(any(array.weights < 0)) {
        stop("array.weights cannot have any negative element")
    }
    if (is.null(array.weights)) {
        array.weights <- -9.9
    } else {
        array.weights <- array.weights / sum(array.weights)
    }

    ## less typing
    first.array <- array.names[1]
    objChrom <- inherits(obj[[first.array]], "RJaCGH.Chrom") 
    objGenome <- inherits(obj[[first.array]], "RJaCGH.Genome") 
    ## we could process objects without chrom as objGenome,
    ## but for coherence with former approaches, and to avoid
    ## ugly "dummyChr" we continue with the objNone approach
    objNone <- ifelse(idchr[1] == dummyChr, TRUE, FALSE)
    if(objNone) objGenome <- FALSE
    
    #####################################################
    ## Main loop: each chromosome is done separately.
    ## For each chrom, do all arrays.
    #####################################################
    
    ## If Genome or None type objects probs of states common for all
    ## chroms.
    if(!objChrom) {
        stretchedProbsList <- lapply(obj[1:narrays],
                                     getStretchedStateProbs)
    }
            
    tmpres <- list()
    for(chromNum in 1:nchrom) {
        if(!objNone) {
            thisChromChar <- idchr[chromNum]
            tmpres[[thisChromChar]] <- list()
        }
        if(verbose) {
            cat("\n -----------------------------------------------")
            cat("\n  Doing chromosome ", chromNum, "\n")
        }
        if(objChrom) {
            nprobes <- obj[[1]][[chromNum]]$lengthY
        } else if(objGenome) {
            nprobes <- sum(obj$Chrom == idchr[chromNum])
        } else {
            nprobes <- length(obj$Chrom) ## FIXME!
            ## should be same as obj[[1]]$lengthY
        }

        if(objChrom) {
            stretchedProbsList <-
                lapply(obj[array.names],
                       function(x) getStretchedStateProbs(x[[chromNum]]))
        }
        filename.tmp <- sapply(obj[array.names],
                               function(x) f2tmp(x$viterbi[[chromNum]]$gzipped_filename,
                                                 tempdir = obj$temp.dir))
        fe <- sapply(filename.tmp, file.exists)
        if(!all(fe)) {
            m1 <- "Some of the gzipped files which should exist don't.\n"
            m2 <- "Please rerun with option 'force.write.files = TRUE'"
            stop(paste(m1, m2))
        }
        filename <- paste(filename.tmp, collapse = "\n")
        
        num.sequences <- sapply(obj[1:narrays],
                                function(x) x$viterbi[[chromNum]]$num_sequences)
        ## Indices as C indices: start at 0.
        ## These give: c(0, last index)
        ## and last index = index of starting pos of a new seq.
        ## so last index = total number of sequences.
        starting.indices.sequences <- c(0, cumsum(num.sequences))
        starting.indices.state.probs <-
            c(0, cumsum(sapply(stretchedProbsList, length)))

        res <- .C("wrap_pREC",
                  alteration = as.integer(alteration.int),
                  numarrays = as.integer(narrays),
                  num_sequences = as.integer(num.sequences),
                  num_probes = as.integer(nprobes),
                  starting_indices_sequences =
                  as.integer(starting.indices.sequences),
                  starting_indices_state_probs =
                  as.integer(starting.indices.state.probs),
                  filename = as.character(filename),
                  threshold = as.double(p),
                  freq_arrays = as.integer(freq.array),
                  array_weights = as.double(array.weights),
                  state_probs = as.double(unlist(stretchedProbsList)),
                  numregions = as.integer(-9),
                  total_narrays = as.integer(-9),
                  regionsStart = as.integer(rep(-9, nprobes)),
                  regionsEnd = as.integer(rep(-9, nprobes)),
                  regionsProb = as.double(rep(-9.9, nprobes)),
                  verboseC = as.integer(ifelse(verbose, 1, 0)),
                  method_prec = as.integer(method_prec))

        if(method == "pREC_S") {
            ## Second call to C, to get results
            number.regionsS <- res$numregions
            total.arrays.S <- res$total_narrays
            rm(res);
            gc()
            if(number.regionsS > 0) {
                res <- .C("return_pREC_S",
                      regionsStart = as.integer(rep(-99, number.regionsS)),
                          regionsEnd = as.integer(rep(-99, number.regionsS)),
                          regionsNarrays = as.integer(rep(-99, number.regionsS)),
                          allArrays = as.integer(rep(-99, total.arrays.S)))
                ## Recall that C gives last elements first, because the head of
                ## of the linked list is the last region found. Thus,
                ## revert order, for easier visualization of results.
                res$regionsStart <- rev(res$regionsStart)
                res$regionsEnd <- rev(res$regionsEnd)
                res$regionsNarrays <- rev(res$regionsNarrays)
                res$allArrays <- rev(res$allArrays)
                res$regionsCumNarrays <- c(0, cumsum(res$regionsNarrays))
                res$numregions <- number.regionsS
                ## Some redundancy above, but leave for now.
                ## could change C code and return cumulative sum directly
                ## as already available from regS->sum_num_arrays.
                ## But messier with reverse order
            } else {
                res <- list()
                res$numregions <- 0
            }
        }
        
        ################ Build output objects   ############
        ### Common to all pREC: get Start and End
        if(!objNone) {
            if (!is.null(obj$Start)) {
                Start <- obj$Start[obj$Chrom == thisChromChar]
                End <- obj$End[obj$Chrom == thisChromChar]
            } else if (!is.null(obj$Pos.rel)) {
                Start <- obj$Pos.rel[obj$Chrom == thisChromChar]
                End <- obj$Pos.rel[obj$Chrom == thisChromChar]
            } else {
                Start <- obj$Pos[obj$Chrom == thisChromChar]
                End <- obj$Pos[obj$Chrom == thisChromChar]
            }
        }

##         if(objChrom) {
##             if (!is.null(obj[[1]][[chromNum]]$Start)) {
##                 Start <- obj[[1]][[chromNum]]$Start
##                 End <- obj[[1]][[chromNum]]$End
##             } else if (!is.null(obj[[1]][[chromNum]]$Pos.rel)) {
##                 Start <- obj[[1]][[chromNum]]$Pos.rel
##                 End <- obj[[1]][[chromNum]]$Pos.rel
##             } else {
##                 Start <- obj[[1]][[chromNum]]$Pos
##                 End <- obj[[1]][[chromNum]]$Pos
##             }
##         } else if(objGenome) {
##             if (!is.null(obj[[1]]$Start)) {
##                 Start <- obj[[1]]$Start[obj[[1]]$Chrom == thisChromChar]
##                 End <- obj[[1]]$End[obj[[1]]$Chrom == thisChromChar]
##             } else if (!is.null(obj[[1]]$Pos.rel)) {
##                 Start <- obj[[1]]$Pos.rel[obj[[1]]$Chrom == thisChromChar]
##                 End <- obj[[1]]$Pos.rel[obj[[1]]$Chrom == thisChromChar]
##             } else {
##                 Start <- obj[[1]]$Pos[obj[[1]]$Chrom == thisChromChar]
##                 End <- obj[[1]]$Pos[obj[[1]]$Chrom == thisChromChar]
##             }
##         }
        if(!objNone) {
            regionstmp <- createRegionsList(res, Start, End,
                                            obj$array.names)
            if(!is.null(regionstmp)) {
                attr(regionstmp, "alteration") <- alteration
                attr(regionstmp, "Chrom") <- thisChromChar
                class(regionstmp) <- paste(method, ".chr", sep = "")
            }
            tmpres[[thisChromChar]] <- regionstmp
        }
    }
    
    if(objNone) { ## outside the chromosome loop as this has no chrom.
        if (!is.null(obj$Start)) {
            Start <- obj$Start
            End <- obj$End
        } else if (!is.null(obj$Pos.rel)) {
            Start <- obj$Pos.rel
            End <- obj$Pos.rel
        } else {
            Start <- obj$Pos
            End <- obj$Pos
        }
        regions <- createRegionsList(res, Start, End, obj$array.names)
        if(is.null(regions)) regions <- list()
        attr(regions, "alteration") <- alteration
        class(regions) <- c(paste(method, ".none", sep = ""), method)
    } else {
        regions <- tmpres
        class(regions) <- c(paste(method, ".Chromosomes", sep = ""),
                            method)
        attr(regions, "alteration") <- alteration
    }
    if(method == "pREC_S") {
        attr(regions, "p") <- p
        attr(regions, "freq.array") <- freq.array
    }
    attr(regions, "array.names") <- array.names
    return(regions)
}    

################ Some auxiliary functions called only from pREC
getStretchedStateProbs <- function(obj) {
    ## for sequence probs.
    ## returns the state probs, from k = 1 to k = max.k,
    ## in row-major order.
    return(unlist(lapply(obj$fit.k[seq_along(levels(obj$k))],
                         function(x) as.vector(t(x$state.labels)))))
}

createRegionsList <- function(res, Start, End, array.names) {
### Creates the "regions" object, a list.
### Used by both pREC-S and pREC-A.
    if(res$numregions == 0)
        return(NULL)
    regions <- list()
    for(i in 1:res$numregions) {
        regions[[i]] <- list()
        regions[[i]]$start <- Start[res$regionsStart[i] + 1]
        regions[[i]]$end <- End[res$regionsEnd[i] + 1]
        regions[[i]]$indexStart <- res$regionsStart[i] + 1
        regions[[i]]$indexEnd <- res$regionsEnd[i] + 1
        regions[[i]]$genes <- (res$regionsEnd[i] - res$regionsStart[i]) + 1
        if(!is.null(res$regionsProb[i])) ##pREC-A
            regions[[i]]$prob <- res$regionsProb[i]
        if(!is.null(res$regionsNarrays[i])) { ## pREC-S
            ## rev for nicer output
            regions[[i]]$arrays <-
                rev(res$allArrays[(res$regionsCumNarrays[i] + 1):
                                  (res$regionsCumNarrays[i + 1])]) + 1
            regions[[i]]$members <- array.names[regions[[i]]$arrays]
        }
    }
    return(regions)
}

f2tmp <- function(x, tempdir) {
    ## the name of tempfiles has "./" under Linux and ".\\" under Windoze.
    ## split, and paste.
    fnam <- paste("rjacgh_seq",
                  unlist(strsplit(x, "rjacgh_seq"))[2], sep = "")
    return(file.path(tempdir, fnam))
}

print.pREC_A.none <- function(x,...) {
    if(length(x) == 0) {
        res <- "No common regions found"
    } else {
        res <- sapply(x, function(y) c(y$start, y$end, y$genes, y$prob)
                      )
        res <- t(res)
        if (ncol(res)>0) {
            colnames(res) <-   c("Start", "End", "Probes",
                                 paste("Prob.", attr(x, "alteration")))
            res <- as.data.frame(res)
        }
        else {
            res <- "No common regions found"
        }
    }
    cat("\n")
    print(res)
}

print.pREC_A <- function(x, ...) {
    UseMethod("print.pREC_A", ...)
}

print.pREC_A.Chromosomes <- function(x,...) {
    if(length(x) == 0) {
        res <- "No common regions found"
    } else {
        res <- NULL
        if(is.null(names(x)))
            names(x) <- seq_along(x)
        for(i in unique(names(x))) {
            if(length(x[[i]]) > 0) {
                tmp <- sapply(x[[i]], function(z) {
                    c(i,  z$start, z$end, z$genes, z$prob)})
                res <- rbind(res, t(tmp))
            }
        }
        if(!is.null(res)) {
            colnames(res) <-   c("Chromosome", "Start", "End", "Probes",
                                 paste("Prob.", attr(x, "alteration")))
            res <- as.data.frame(res)
        }
    }
    cat("\n")
    print(res)
}

print.pREC_S <- function(x, ...) {
    UseMethod("print.pREC_S", ...)
}

print.pREC_S.none <- function(x,...) {
    cat("Common regions of", attr(x, 'alteration'), "of at least",
        attr(x, 'p'), "probability:\n")
    if(length(x) == 0) {
        res <- "No regions found"
        print(res)
    } else {
        res <- sapply(x, function(z)  {
            c(z$start, z$end, z$indexEnd - z$indexStart + 1,
              paste(z$members, collapse=";"))
            
        })
        res <- t(res)
        if (!is.null(res)) {
            colnames(res) <-   c("Start", "End", "Probes",
                                 "Arrays")
            res <- as.data.frame(res)
            print(res)
        }
    }
}


print.pREC_S.Chromosomes <- function(x,...) {
    cat("Common regions of", attr(x, 'alteration'), "of at least",
        attr(x, 'p'), "probability:\n")
    if(length(x) == 0) {
        res <- "No regions found"
        print(res)
        return()
    } else {
        res <- NULL
        if (is.null(names(x)))
            names(x) <- 1:length(x)
        for(i in unique(names(x))) {
            if(length(x[[i]]) > 0) {
                tmp <- sapply(x[[i]], function(z)  {
                    c(i, z$start, z$end, z$indexEnd - z$indexStart + 1,
                      paste(z$members, collapse=";"))
                })
                
                res <- rbind(res, t(tmp))
            }
        }
    }
    if (!is.null(res)) {
        colnames(res) <-   c("Chromosome", "Start", "End", "Probes",
                             "Arrays")
        res <- as.data.frame(res)
        print(res)
    }
}

plot.pREC_S <- function(x, array.labels=NULL,
                        stats=TRUE,
                        col=NULL, breaks=NULL,
                        dend=TRUE, method="single",
                        Chrom=NULL, ...) {
    #### Helper functions used only here
    f.create.grid <- function(x, array.names) {
        probes <- x$indexStart:x$indexEnd
        members <- factor(x$members, levels=array.names)
        probes.length <- (x$end - x$start + 1) / length(probes)
        expand.grid(members, members, probes, probes.length)
    }
    f2 <- function(z, array.names) {
        tmp <- lapply(z,  function(x)
                      f.create.grid(x, array.names = array.names))
        tmp <- do.call("rbind", tmp)
        tmp <- unique(tmp)
        return(list(internal.inc.mat = table(tmp[ , -c(3,4)]),
                    internal.length.mat =
                    xtabs(Var4 ~ Var1 + Var2, data=tmp)))
    }


    array.names <- attr(x, 'array.names')
    if (is.null(array.labels)) array.labels <- array.names
    k <- length(array.names)
    par(oma=c(2, 2, 4, 2))
    layout(matrix(c(1, 2), 1, 2), widths=c(1, 7))
    if(!is.null(Chrom)) {
        if(!inherits(x, "pREC_S.Chromosomes")) {
            stop("With this type of object it makes no sense to use the Chrom argument")
        }
        ## Don't do this before, as array.names
        ## not available in x[[Chrom]]
        x <- x[[Chrom]]
        if(is.null(x)) {
            stop("No common regions with Chromosome ", Chrom)
        }
        class(x) <- "pREC_S.none"
    }

    ###  The only difference between objects with or without chromosomes:
    ###  this if
    if(inherits(x, "pREC_S.none")) {
        t1 <- f2(x, array.names)
        inc.mat <- t1$internal.inc.mat
        length.mat <- t1$internal.length.mat
    } else {
        inc.mat <- matrix(0, length(array.names), length(array.names))
        length.mat <- matrix(0, length(array.names), length(array.names))

        for(chr in names(x)) {
            if(length(x[[chr]]) > 0) {
                t1 <- f2(x[[chr]], array.names)
                inc.mat <- inc.mat + t1$internal.inc.mat
                length.mat <- length.mat + t1$internal.length.mat
            }
        } 
    }
    #####  The rest is all common
    
    length.mat[inc.mat > 0] <- length.mat[inc.mat > 0] /
        inc.mat[inc.mat > 0]
    diag(inc.mat) <- 0
    diag(length.mat) <- 0
    distances <- 1 - (inc.mat / matrix(max(inc.mat),
                                       nrow(inc.mat),
                                       ncol(inc.mat)))
    obj.dend <- as.dendrogram(hclust(as.dist(distances),
                                     method=method))
    par(mai=c(0.5, 0, 0.5, 0.5))
    if (dend) {
        reordering <- order.dendrogram(obj.dend)
        inc.mat <- inc.mat[reordering, reordering]
        length.mat <- length.mat[reordering, reordering]
        plot(obj.dend, ylab="", main="", axes=FALSE,
             xlab="", horiz=TRUE, yaxs="i", leaflab="none")
    }
    else {
        reordering <- 1:length(array.labels)
        plot.new()
    }
    if (is.null(col)) {
        ## default palette taken from redgreen from
        ## gplots, Gregory R. Warnes
        col <- c("#FF0000", "#DF0000", "#BF0000", "#800000", "#600000",
                 "#400000", "#200000", "#000000",  "#002000", "#004000",
                 "#006000", "#008000", "#00BF00", "#00DF00", "#00FF00")
        if (attr(x, "alteration") == "Gain") {
            col <- col[8:1]
        }
        else {
            col <- col[8:15]
        }
    }
    if (is.null(breaks)) {
        breaks <- quantile(inc.mat[inc.mat>0],
                           p=seq(from=0, to=1,
                           length=length(col) - 1))
        breaks <- c(0, 0.1, breaks)
        breaks <- breaks - 0.05
        breaks[length(breaks)] <- breaks[length(breaks)] + 0.10
    }
    
    image(x=1:k, y=1:k, z=inc.mat,
          axes=FALSE, col=col, breaks=breaks, xlab="", ylab="",...)
    axis(side=1, at=1:k, labels=array.labels[reordering], las=2,
         tick=FALSE,...)
    axis(side=2, at=1:k, labels=array.labels[reordering],
         las=2, tick=FALSE,...)
    box()
    if (stats) {
        for(i in 1:k) {
            text(rep(i, k), 1:k,
                 paste(inc.mat[,i], " (", round(length.mat[,i], 1), ")",
                       sep=""), col="white", ...)
        }
    }
    list(probes=inc.mat, length=length.mat)
}


getHostname.System <- function (static, ...) {
    ## From the function of the same name in package R.utils (v. 1.0.1)
    ## by Henrik Bengtsson
    host <- Sys.getenv(c("HOST", "HOSTNAME", "COMPUTERNAME"))
    host <- host[host != ""]
    if (length(host) == 0) {
        host <- Sys.info()["nodename"]
        host <- host[host != ""]
        if (length(host) == 0) {
            host <- readLines(pipe("/usr/bin/env uname -n"))
        }
    }
    host[1]
}

the.time.with.ms <- function() {
    uu <- as.POSIXlt(Sys.time())
    return(paste(uu$hour, uu$min,
                 paste(unlist(strsplit(as.character(uu$sec), "\\.")),
                       collapse = ""), sep = ""))
}

tempdir2 <- function() {
    direxists <- TRUE
    while(direxists) {
        p1 <-  paste(getHostname.System(), round(runif(1, 1, 9999)),
                     the.time.with.ms(), sep = "_")
        p1 <- paste(tempfile(pattern = "tmpdir_RJaCGH_",
                             tmpdir = "."),
                    p1, sep = "_")
        if(!file.exists(p1)) direxists <- FALSE
    }
    dir.create(p1)
    return(p1)
}

tempfile2 <- function(pattern = "rjacgh_seq", tmpdir = ".") {
    ## return complete path name to a temp file
    ## generates the root filename
    ## no file type accepted
    p00 <- tempfile00(pattern = pattern)
    filename <- tempfile(pattern = p00, tmpdir = tmpdir)
    return(filename)
}

tempfile3 <- function(fname, tmpdir = ".", ftype = "RData") {
    ## return complete path name to a temp file
    ## root filename is parameter
    filename <- paste(
                      tempfile(pattern = fname, tmpdir = tmpdir),
                      ftype, sep = ".")
    return(filename)
}

tempfile00 <- function(pattern = "rjacgh_seq") {
    ## generates "root" (no dir) of a temp file
    pattern1 <- paste(getHostname.System(), round(runif(1, 1, 9999)),
                      the.time.with.ms(), sep = "_")
    return(paste(pattern, pattern1, sep = "_"))
}





### Using save and load or saveRDS and readRDS are about equally
### fast and memory consuming. But the RDS files deal with an object
### directly, not the name of an object. Some experiments follow way down.

writeRO2 <- function(x, rootn, tmpdir, ftype = "RData") {
    fname <- tempfile3(tempfile00(rootn), ".", ftype = ftype)
    fnfull <- file.path(tmpdir, fname)
    saveRDS(x, file = fnfull, compress = FALSE)
    return(list(tmpdir = tmpdir, fname = fname))
}

readRO2 <- function(x) {
    fname <- file.path(x$tmpdir, x$fname)
    if(!file.exists(fname)) {
        stop(paste("Cannot acces file", fname,
                   ". Either you are in the wrond directory",
                   "or the file has been deleted (the later is ",
                   "a non-recoverable error)."))
    }
    return(readRDS(fname))
}


### writeRO <- function(x, rootn, tmpdir) {
###     ## x has to be an object name
### ############################################
###   #########################################
###   ## tempfile2 has only 2 arguments
###     fn <- tempfile2(rootn, tmpdir, ftype = "RData")
###     save(x, file = fn, compress = FALSE)
###     ## zz: FIXME: broken! does not return filenames
### }

### readRO <- function(x) {
###     fname <- file.path(x$tmpdir, x$fname)
###     if(!file.exists(fname)) {
###         stop(paste("Cannot acces file", fname,
###                    ". Either you are in the wrond directory",
###                    "or the file has been deleted (the later is ",
###                    "a non-recoverable error)."))
###     }
###     return(get(load(fname)))
### }





###################################################################
###################################################################

######       chainsSelect stuff     ###############################




###  I stop work with that. The problem is how to, later, collapse
###  the Viterbi sequences.

### Guide to the new code:
### try to have a single select chains (chainsSelect.New)
### that iterates slightly differntly depending on model Chrom or genome
### two small new functs. are used


### gelman.rubin is also obsoleted (does not apply) and collapseChains too.

### note some comments on "multipleChains" function above for further work.



### chainsSelect.New <- function(obj, nutrim = NULL, trim = NULL) {
###     if((is.null(nutrim) & is.null(trim)) |
###        (!is.null(nutrim) & !is.null(trim)))
###         stop("Exactly one of nutrim or trim must have non-NULL values")
    
###     n <- length(obj) ## obj is a list
###     uChrom <- unique(obj[[1]]$Chrom)
###     array.names <- obj[[1]]$array.names
###     if (is.null(nutrim)) {
###         nutrim <- floor(n * trim)
###         trim <- NULL
###     }
###     LH <- selectChainsLimit(n, nutrim, trim)
###     newobj <- list()
###     for (j in 1:(n - nutrim)) {
###         newobj[[j]] <- list()
###         for(arr in array.names) 
###             newobj[[j]][[arr]] <- list()
###     }
###     viterbi.tmp <- list()
###     for (arr in array.names) {
###         viterbi.tmp[[arr]] <- list()
###         if(model == "Chrom") {
###             for (chr in uChrom) {
###                 viterbi.tmp[[arr]][[chr]] <- list()
###                 cssl <-
###                     chainsToRetain(lapply(obj[1:n],
###                                           function(x) x[[arr]][[chr]]$k),
###                                    LH['Low'], LH['High'])
###                 for(j in 1:(n - nutrim)) {
###                     newobj[[j]][[arr]][[chr]] <- obj[[ cssl[j] ]][[arr]][[chr]]
###                     viterbi.tmp[[arr]][[chr]][[j]] <-
###                         obj[[ cssl[j] ]][[arr]]$viterbi[[chr]]
###                 }
###             }
###             for(j in 1:(n - nutrim)) {
###                 for(chr in uChrom) {
###                     newobj[[j]][[arr]]$viterbi[[chr]] <-
###                         viterbi.tmp[[arr]][[chr]][[j]]
###                 }
###             }
###         } else { ##model is Genome
###             cssl <- chainsToRetain(lapply(obj[1:n], function(x) x[[arr]]$k),
###                                    LH['Low'], LH['High'])
###             for(j in 1:(n - nutrim)) {
###                 newobj[[j]][[arr]] <- obj[[ cssl[j] ]][[arr]]
###             }
###         }
###     }
            
            
### ###         obj.temp <- list()
### ### ###         for (j in 1:n) {
### ### ###             obj.temp[[j]] <- obj[[j]][[i]]
### ### ###         }

### ###         obj.temp <- chainsSelect(obj.temp, nutrim = nutrim, trim = trim,
### ###                                  uChrom = uChrom)
### ###         for (j in 1:(n-nutrim)) {
### ###             newobj[[j]][[i]] <- obj.temp[[j]]
### ###             class(newobj[[j]]) <- class(obj[[j]])
### ###         }
### ###     }
###     for (j in 1:(n-nutrim)) {
###         newobj[[j]]$array.names <- obj[[j]]$array.names
###         newobj[[j]]$Pos <- obj[[j]]$Pos
###         newobj[[j]]$Chrom <- obj[[j]]$Chrom
###         newobj[[j]]$Dist.for.model <- obj[[j]]$Dist.for.model
###         newobj[[j]]$call <- obj[[j]]$call
###         newobj[[j]]$temp.dir <- obj[[j]]$temp.dir
###         ##    attr(newobj[[j]], 'names') <- attr(obj[[j]], 'names')
###     }
###     class(newobj) <- class(obj)
###     newobj
### }


### selectChainsLimits <- function(n, nutrim, trim) {
###     if((is.null(nutrim) & is.null(trim)) |
###        (!is.null(nutrim) & !is.null(trim)))
###         stop("Exactly one of nutrim or trim must have non-NULL values")

###     if (is.null(nutrim)) {
###         nutrim <- round(trim * n)
###         if(nutrim >= n)
###             stop(paste("Specify a smaller trim: your trim ",
###                        trim, "leads to selecting no chains"))
###     }
###     else {
###         if(nutrim >= n)
###             stop("Number to trim must be smaller than number of chains")
###     }
###     lo <- floor(nutrim/2) + 1
###     hi <- n - (nutrim - lo + 1)
###     return(c(Low = lo, High = hi))
### }


### chainsToRetain <- function(obj, lo, hi) {
###     ##obj: a list; each entry is a the "k"
###     cat("\n chains to retain \n")
###     meank <- unlist(lapply(obj,function(x)
###                            mean(as.numeric(as.character(x)))))
###     keepInd <- order(meank)[lo:hi]
###     print(paste("keepInd ", paste(keepInd, collapse = " ")))
###     return(keepInd)
### }





### chainsSelect <- function(obj, nutrim = NULL, trim = NULL, ...) {
###   UseMethod("chainsSelect", obj[[1]])
### }

### chainsSelect.iRJaCGH <- function(obj, nutrim = NULL, trim = NULL, ...) {
###     cat("\ncs iRJ\n")

###     ## This uses too much memory: obj is passed by
###     ## copy, and we make yet another partial copy at the end
###     n <- length(obj)
###     if((is.null(nutrim) & is.null(trim)) |
###        (!is.null(nutrim) & !is.null(trim)))
###         stop("Exactly one of nutrim or trim must have non-NULL values")

###     if (is.null(nutrim)) {
###         nutrim <- round(trim * n)
###         if(nutrim >= n)
###             stop(paste("Specify a smaller trim: your trim ",
###                        trim, "leads to selecting no chains"))
###     }
###     else {
###         if(nutrim >= n)
###             stop("Number to trim must be smaller than number of chains")
###     }
###     lo <- floor(nutrim/2) + 1
###     hi <- n - (nutrim - lo + 1)
###     meank <- unlist(lapply(obj,function(x)
###                            mean(as.numeric(as.character(x$k)))))
###     keepInd <- order(meank)[lo:hi]
###     print(paste("keepInd ", paste(keepInd, collapse = " ")))
###     newobj <- list()
###     j <- 1
###     for(i in keepInd) {
###         newobj[[j]] <- obj[[i]]
###         j <- j + 1
###     }
###     class(newobj) <- class(obj)
###     newobj
### }





### chainsSelect.RJaCGH.Genome <- function(obj, nutrim = NULL, trim = NULL, uChrom) {
###     cat("\ncs genome\n")
###     newobj <- chainsSelect.iRJaCGH(obj, nutrim = nutrim, trim = trim, uChrom)
###     class(newobj) <- class(obj)
###     newobj
### }

### chainsSelect.RJaCGH.Chrom <- function(obj, nutrim = NULL, trim = NULL, uChrom) {
###     cat("\ncs chrom\n")
###     if((is.null(nutrim) & is.null(trim)) |
###        (!is.null(nutrim) & !is.null(trim)))
###         stop("Exactly one of nutrim or trim must have non-NULL values")
###     newobj <- list()
###     n <- length(obj)
###     if (is.null(nutrim)) {
###         nutrim <- floor(n * trim)
###         trim <- NULL
###     }
###     ## trim
###     for (i in 1:(n-nutrim)) {
###         newobj[[i]] <- list()
###     }
###     for (chr in uChrom) {
###         oldobj <- lapply(obj, function(x) x[[chr]])
###         tmpobj <- chainsSelect.iRJaCGH(oldobj, nutrim=nutrim, trim=trim)
###         for (i in 1:(n-nutrim)) {
###             newobj[[i]][[chr]] <- tmpobj[[i]]
###         }
###     }
###     for (i in 1:(n-nutrim)) {
###         newobj[[i]]$viterbi <- obj[[i]]$viterbi
###         ##      newobj[[i]]$Pos <- obj[[i]]$Pos
###         ##      newobj[[i]]$Pos.rel <- obj[[i]]$Pos.rel
###         ##      newobj[[i]]$Chrom <- obj[[i]]$Chrom
###         ## que hacia esto? zz: FIXME
###         ##     attr(newobj[[i]], 'names') <- attr(obj[[i]], 'names')
###         class(newobj[[i]]) <- "RJaCGH.Chrom"
###     }
###     newobj
### }



### chainsSelect.formerdefault <- function(obj, nutrim = NULL, trim = NULL) {
###   if((is.null(nutrim) & is.null(trim)) |
###      (!is.null(nutrim) & !is.null(trim)))
###     stop("Exactly one of nutrim or trim must have non-NULL values")
  
###   n <- length(obj) ## obj is a list
###   newobj <- list()
###   uChrom <- unique(obj[[1]]$Chrom)
###   array.names <- obj[[1]]$array.names
###   if (is.null(nutrim)) {
###     nutrim <- floor(n * trim)
###     trim <- NULL
###   }
###   for (j in 1:(n-nutrim)) {
###     newobj[[j]] <- list()
###   }
###   for (i in array.names) {
###     obj.temp <- list()
###     for (j in 1:n) {
###       obj.temp[[j]] <- obj[[j]][[i]]
###     }
###     obj.temp <- chainsSelect(obj.temp, nutrim = nutrim, trim = trim,
###                             uChrom = uChrom)
###     for (j in 1:(n-nutrim)) {
###       newobj[[j]][[i]] <- obj.temp[[j]]
###       class(newobj[[j]]) <- class(obj[[j]])
###     }
###   }
###   for (j in 1:(n-nutrim)) {
###     newobj[[j]]$array.names <- obj[[j]]$array.names
###     newobj[[j]]$Pos <- obj[[j]]$Pos
###     newobj[[j]]$Chrom <- obj[[j]]$Chrom
###     newobj[[j]]$Dist.for.model <- obj[[j]]$Dist.for.model
###     newobj[[j]]$call <- obj[[j]]$call
###     newobj[[j]]$temp.dir <- obj[[j]]$temp.dir
### ##    attr(newobj[[j]], 'names') <- attr(obj[[j]], 'names')
###   }
###   class(newobj) <- class(obj)
###   newobj
### }




### ## zz: FIXME: RETHING HOW WE STORE MULTIPLE CHAINS!!!
### ## very wasteful in terms of storage to use a list of RJaCGH obj.
### ## and slow access
### ## Now a bit less, but still
### ## try also to support futhre parallelization

### multipleChains <- function(...) {
###     ## "standard" RJaCGH
###     ## the rest as lists, with only unique stuff.
###     ## change how RJaCGH runs

###     ## New structure ought to ease:
###     ##  a) collapseChain
###     ##  b) chainsSelect [change to selectChain]
###     ##  c) gelman.rubin plot
    
###     ## Common:
###     ## "array.names"    "Pos"            "Chrom"         
###     ## "Dist.for.model" "call"           "temp.dir"      

###     ## other things not to keep: Y.
    
###     ## BEWARE: temp.dir!!!
### }
    


### collapseChain <- function(obj) {
###   UseMethod("collapseChain", obj[[1]])
### }

### collapseChain.iRJaCGH <- function(obj) {
###   normal.reference <- attr(obj, "normal.reference")
###   window <- attr(obj, "window")
###   singleState <- attr(obj, "singleState")
    
###   newobj <- list()
###   class(newobj) <- "RJaCGH"
###   newobj$y <- NULL
###   newobj$x <- NULL
###   newobj$Pos <- NULL
###   newobj$prob.b <- NULL
###   newobj$prob.d <- NULL
###   newobj$prob.s <- NULL
###   newobj$prob.c <- NULL
###   newobj$k <- NULL
###   C <- length(obj)
###   k <- max(as.numeric(as.character(levels(obj[[1]]$k))))
###   for (i in 1:k) {
###     newobj$fit.k[[i]] <- list()
###     newobj$fit.k[[i]]$stat <- NULL
###     newobj$fit.k[[i]]$mu <- NULL
###     newobj$fit.k[[i]]$sigma.2 <- NULL
###     newobj$fit.k[[i]]$beta <- NULL
###     newobj$fit.k[[i]]$loglik <- NULL
###     newobj$fit.k[[i]]$prob.states <- matrix(0, nrow=length(obj[[1]]$y),
###   ncol=i)
### ###     newobj$fit.k[[i]]$state.labels <- NULL
###   }
###   newobj$prob.b <- c(0, 0)
###   newobj$prob.d <- c(0, 0)
###   newobj$prob.s <- 0
###   newobj$prob.c <- 0
###   for (i in 1:C) {
###     newobj$k <- c(newobj$k, obj[[i]]$k)
###     newobj$prob.b <- newobj$prob.b + obj[[i]]$prob.b
###     newobj$prob.d <- newobj$prob.d + obj[[i]]$prob.d
###     newobj$prob.s <- newobj$prob.s + obj[[i]]$prob.s
###     newobj$prob.c <- newobj$prob.c + obj[[i]]$prob.c
###     newobj$y <- obj[[i]]$y
###     newobj$x <- obj[[i]]$x
###     newobj$Pos <- obj[[i]]$Pos
###     newobj$Pos.rel <- obj[[i]]$Pos.rel
###     for (j in 1:k) {
###       newobj$fit.k[[j]]$stat <- obj$fit.k[[i]][[j]]$stat
###       newobj$fit.k[[j]]$mu <- rbind(newobj$fit.k[[j]]$mu, obj$fit.k[[i]][[j]]$mu)
###       newobj$fit.k[[j]]$sigma.2 <- rbind(newobj$fit.k[[j]]$sigma.2, obj$fit.k[[i]][[j]]$sigma.2)
###       ## FIXME: do not grow things this way!!!
###       newobj$fit.k[[j]]$beta <- c(newobj$fit.k[[j]]$beta, obj$fit.k[[i]][[j]]$beta)
###       newobj$fit.k[[j]]$loglik <- c(newobj$fit.k[[j]]$loglik,
###                               obj$fit.k[[i]][[j]]$loglik)

     
###       if (!is.null(obj$fit.k[[i]][[j]]$prob.states)) {
###         newobj$fit.k[[j]]$prob.states <- newobj$fit.k[[j]]$prob.states + 
###           obj$fit.k[[i]][[j]]$prob.states * nrow(obj$fit.k[[i]][[j]]$mu)
###       }
###       else {
###         newobj$fit.k[[j]]$prob.states <- newobj$fit.k[[j]]$prob.states
###       }
###     }
###   }
###   for (j in 1:k) {
###     newobj$fit.k[[j]]$beta <- array(newobj$fit.k[[j]]$beta,
###                               c(j, j, length(newobj$fit.k[[j]]$beta)/(j*j)))
###     newobj$fit.k[[j]]$prob.mu <- length(unique(newobj$fit.k[[j]]$mu[,1])) /
###       length(newobj$fit.k[[j]]$mu[,1])
###     newobj$fit.k[[j]]$prob.sigma.2 <- length(unique(newobj$fit.k[[j]]$sigma.2[,1])) /
###       length(newobj$fit.k[[j]]$sigma.2[,1])
###     if (j >1) {
###       newobj$fit.k[[j]]$prob.beta <- length(unique(newobj$fit.k[[j]]$beta[1,2,])) /
###         length(newobj$fit.k[[j]]$beta[1,2,])
###     }
###     else newobj$fit.k[[j]]$prob.beta <- NULL
###     if(nrow(newobj$fit.k[[j]]$mu)) {
###       newobj$fit.k[[j]]$prob.states <- newobj$fit.k[[j]]$prob.states /
###         nrow(newobj$fit.k[[j]]$mu)
###     }
###     else {
###       newobj$fit.k[[j]]$prob.states <- NULL
###     }
###   }
###   newobj$k <- factor(newobj$k, levels=1:k)
###   ## Recompute state labels
###   sllist <- relabel.core(obj = newobj,
###                          normal.reference = normal.reference,
###                          window = window,
###                          singleState = singleState)
###   k.max <- max(as.numeric(levels(newobj$k)))
###   for(i in (1:k.max))
###       newobj$fit.k[[i]]$state.labels <- sllist[[i]]
 
###   attr(newobj, "normal.reference") <- normal.reference
###   attr(newobj, "window") <- window
###   attr(newobj, "singleState") <- singleState

###   newobj
### }



### collapseChain.RJaCGH.Genome <- function(obj) {
###   newobj <- collapseChain.RJaCGH(obj)
###   newobj$model <- obj[[1]]$model
###   newobj$Chrom <- obj[[1]]$Chrom
###   class(newobj) <- class(obj[[1]]) ## zz: fixme: will this work?
###   newobj
### }

### collapseChain.RJaCGH.Chrom <- function(obj) {
###   C <- length(obj)
###   newobj <- list()
###   for (i in unique(obj[[1]]$Chrom)) {
###     obj.temp <- list()
###     for (j in 1:C) {
###       obj.temp[[j]] <- obj[[j]][[i]]
###     }
###     newobj[[i]] <- collapseChain(obj.temp)
###   }
###   newobj$model <- obj[[1]]$model
###   newobj$Chrom <- obj[[1]]$Chrom
###   newobj$Pos <- obj[[1]]$Pos
###   newobj$Pos.rel <- obj[[1]]$Pos.rel
###   class(newobj) <- "RJaCGH.Chrom"
###   newobj
### }

### collapseChain.RJaCGH <- function(obj) {
###   C <- length(obj)
###   newobj <- list()
###   newobj$array.names <- obj[[1]]$array.names
###   for (i in newobj$array.names) {
###     newobj[[i]] <- list()
###     obj.temp <- list()
###     for (j in 1:C) {
###       obj.temp[[j]] <- obj[[j]][[i]]
###     }
###     newobj[[i]] <- collapseChain(obj.temp)
###   }
  
###   class(newobj) <- "RJaCGH"
###   newobj
### }



### ##Gelman, A. and Rubin, D.
### ##Inference from Iterative simulation using multiple sequences
### ##Statistical Science 7, 457--511

### gelman.rubin.plot <- function(obj, bin=1000, array=NULL, Chrom=NULL, k=NULL) {
###   obj.R <- list()
###   if (class(obj[[1]])=="RJaCGH.array" && is.null(array)) {
###     stop("Must specify array\n")
###   }
###   if (class(obj[[1]])=="RJaCGH.Chrom" && is.null(Chrom)) {
###     stop("Must specify chromosome\n")
###   }
  
###   for (i in 1:length(obj)) {
###     if (!is.null(array)) obj[[i]] <- obj[[i]][[array]]
###     if (!is.null(Chrom)) obj[[i]] <- obj[[i]][[Chrom]]
###   }
###   par(mfrow=c(2,2), oma=c(1,1,3,1))

###   C <- length(obj)
###   ## Graph of number of states
###   T <- min(unlist(lapply(obj, function(x) length(x$k))))
###   t <- floor(T / bin)
###   R <- rep(0, t+1)
  
###   for (i in 1:(t+1)) {

###     chain.mean <- rep(0, C)
###     chain.var <- rep(0, C)
###     for(j in 1:C) {
###       if (i>t) batch <- obj[[j]]$k
###       else batch <- obj[[j]]$k[1:(bin*i)]
###       chain.mean[j] <- mean(as.numeric(as.character(batch)))
###       chain.var[j] <- var(as.numeric(as.character(batch)))    
###     }
###     n <- length(batch)
###     B.k <- sum(((chain.mean-mean(chain.mean))^2) * (n)/(C-1))
###     W.k <- mean(chain.var)
###     R[i] <- sqrt((((n-1) / n) * W.k + B.k/n) / W.k)
###   }
###   obj.R$k <- R[length(R)]
###   plot(seq(bin, by=bin, length=t+1), R, type="l",
###        ylim=c(0.8, 1.2), main="Parameter: k", xlab="Batch size")
###   abline(h=1, lty=2)
###   abline(h=c(0.9, 1.1), lty=3)

###   ## Choose k, the max of all chains (changed 07/07/08)
###   if (is.null(k)) {
###     k <- apply(do.call(rbind, lapply(obj, function(x) summary(x$k))), 2, sum)
###     k <- as.numeric(names(which.max(k)))
###   }
###   for (i in 1:C)   obj[[i]] <- obj[[i]][[k]]
  
###   T <- min(unlist(lapply(obj, function(x) nrow(x$mu))))
###    if (T==0)
###      stop ("MCMC not converged.\n Gelman-Brooks values can't be computed and graph will not be drawn")
###   ## no two chains ever visited the same state
###   t <- floor(T / bin)

###   ## Graph of mu  
###   R <- matrix(0, t+1, k)
###   for (i in 1:(t+1)) {
###     chain.mean <- matrix(0, C, k)
###     chain.var <- matrix(0, C, k)
###     ## Number of observations in every chain
###     TOT <- rep(0, C)
###     for (j in 1:C) {
###       if (i>t) batch <- obj[[j]]$mu
###       else batch <- obj[[j]]$mu[1:(bin*i), ]
###       batch <- matrix(batch, ncol=k)
###       TOT[j] <- nrow(batch)
###       chain.mean[j,] <- colMeans(batch)
###       chain.var[j,] <- apply(batch, 2, var)
###     }

###     mean.chain.mean <- matrix(colMeans(chain.mean), C, k, byrow=TRUE)
###     B <- ((chain.mean - mean.chain.mean)^2) * matrix(TOT/(C-1), C, k)
###     B <- colSums(B)
###     W <- colMeans(chain.var)
###     n <- mean(TOT)
###     R[i,] <- sqrt((((n-1) / n) * W + B/n) / W)
###   }
###   matplot(seq(bin, by=bin, length=t+1), R, type="l",
###           ylim=c(0.8, 1.2), main="Parameter: mean", xlab="Batch size")
###   abline(h=c(0.9, 1.1), lty=3)
###   abline(h=1, lty=2)
###   obj.R$mu <- R[nrow(R),]
###   ## Graph of sigma2
###   R <- matrix(0, t+1, k)
###   for (i in 1:(t+1)) {
###     chain.mean <- matrix(0, C, k)
###     chain.var <- matrix(0, C, k)
###     ## Number of observations in every chain
###     TOT <- rep(0, C)
###     for (j in 1:C) {
###       if (i>t) batch <- obj[[j]]$sigma.2
###       else batch <- obj[[j]]$sigma.2[1:(bin*i), ]
###       batch <- matrix(batch, ncol=k)
###       TOT[j] <- nrow(batch)
###       chain.mean[j,] <- colMeans(batch)
###       chain.var[j,] <- apply(batch, 2, var)
###     }
###     mean.chain.mean <- matrix(apply(chain.mean, 2, mean), C, k, byrow=TRUE)
###     B <- ((chain.mean - mean.chain.mean)^2) * matrix(TOT/(C-1), C, k)
###     B <- colSums(B)
###     W <- colMeans(chain.var)
###     n <- mean(TOT)
###     R[i,] <- sqrt((((n-1) / n) * W + B/n) / W)

###   }
###   matplot(seq(bin, by=bin, length=t+1), R, type="l",
###           ylim=c(0.8, 1.2),main="Parameter: sigma.2", xlab="Batch size")
###   abline(h=1, lty=2)
###   abline(h=c(0.9, 1.1), lty=3)
###   obj.R$sigma.2 <- R[nrow(R),]
###   ## Graph of beta
###   if (k >1) {
###     R <- matrix(0, t+1, k*(k-1))  
###     for (i in 1:(t+1)) {
###       chain.mean <- matrix(0, C, k*(k-1))
###       chain.var <- matrix(0, C, k*(k-1))

###       ## Number of observations in every chain
###       TOT <- rep(0, C)
###       for (j in 1:C) {
###         batch <-NULL
###         if (i>t) {
###           for (r in 1:k) {
###             for (c in 1:k) {
###               if (r!=c) batch <- cbind(batch,log(obj[[j]]$beta[r,c,]))
###             }
###           }
###         }
###         else {
###           for (r in 1:k) {
###             for (c in 1:k) {
###               if (r!=c) batch <- cbind(batch,log(obj[[j]]$beta[r, c, 1:(bin*i)]))
###             }
###           }
###         }
###         TOT[j] <- nrow(batch)
###         chain.mean[j,] <- colMeans(batch)
###         chain.var[j,] <- apply(batch, 2, var)
###       }
###       mean.chain.mean <- matrix(apply(chain.mean, 2, mean), C, k*(k-1), byrow=TRUE)
###       B <- ((chain.mean - mean.chain.mean)^2) * matrix(TOT/(C-1), C, k*(k-1))
###       B <- colSums(B)
###       W <- colMeans(chain.var)
###       n <- mean(TOT)
###       R[i,] <- sqrt((((n-1) / n) * W + B/n) / W)
###     }
###     matplot(seq(bin, by=bin, length=t+1), R, type="l",
###             ylim=c(0.8, 1.2), main="Parameter: beta", xlab="Batch size")
###     abline(h=1, lty=2)
###     abline(h=c(0.9, 1.1), lty=3)
###   }
###   obj.R$beta <- R[nrow(R),]
###   mtext("Gelman-Rubin diagnostic plots", outer=TRUE)
###   obj.R
### }






## Examples for chainsSelect
###### A test of chainsSelect

## one array

## y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##        rnorm(100,0, 1)) 
## Pos <- runif(230)
## Pos <- cumsum(Pos)
## Chrom <- rep(1:23, rep(10, 23))

## jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
##            sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5, tau.split.beta=0.5)

## fit.genome <- list()
## for (i in 1:8) {
##     fit.genome[[i]] <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }


## fit.chrom <- list()
## for (i in 1:8) {
##     fit.chrom[[i]] <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }





## ## many arrays

## y1 <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##        rnorm(100,0, 1)) 
## y2 <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##        rnorm(100,0, 1)) 
## y3 <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##         rnorm(100,0, 1)) 
## Y <- cbind(y1, y2, y3)


## fit.genome.arrays <- list()
## for (i in 1:8) {
##     fit.genome.arrays[[i]] <- RJaCGH(y=Y, Pos=Pos, Chrom=Chrom, model="genome",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }


## fit.chrom.arrays <- list()
## for (i in 1:8) {
##     fit.chrom.arrays[[i]] <- RJaCGH(y=Y, Pos=Pos, Chrom=Chrom, model="Chrom",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }


## #####

## o1 <- chainsSelect(fit.chrom.arrays)
## o2 <- chainsSelect(fit.genome.arrays)
## o3 <- chainsSelect(fit.chrom)
## o4 <- chainsSelect(fit.genome)








##############################################################
##############################################################
##############################################################


## Not correct (does not take into account mixing proportions)

## RJaCGHDensity <- function(obj, k=NULL,...) {
##   UseMethod("RJaCGHDensity")
## }

## RJaCGHDensity.RJaCGH <- function(obj, k=NULL,...) {
##   hist(obj$y, probability=TRUE,...)
##   if (is.null(k)) {
##     k <- which.max(summary(obj$k))
##   }
##   obj.sum <- summary(obj, k=k, quantiles=0.5)
##   x <- seq(from=min(obj$y), to=max(obj$y), length=1000)
##   y <- matrix(0, 1000, k)
##   for (i in 1:k) {
##     y[,i] <- 
##       dnorm(x, obj.sum$mu[i], sqrt(obj.sum$sigma.2[i]))
##   }
##   y <- apply(y, 1, mean)
##   lines(x, y, col=2)
## }








####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################


######### Old code, no longer used



### getSequence <- function(obj, filename, alteration) {

###   K <- max(as.numeric(as.character(levels(obj$k))))
###   probs <- prop.table(table(obj$k))
  
###   for (k in 1:K) {
###     if(probs[k] > 0) {
###       inds <- grep(paste("[", substr(alteration, 1, 1), "]"),
###                    obj[[k]]$state.labels)
###       if (length(inds) > 0) {
###         ## Viterbi of all iterations
###         N <- nrow(obj[[k]]$mu)
###         if (inherits(obj, "RJaCGH.genome")) {
###           index <- c(which(!duplicated(obj$Chrom)) - 1, length(obj$y))
###           genome <- length(index) -1
###         }
###         else {
###           index <- c(0, length(obj$y))
###           genome <- 1
###         }
###         obj$x[is.na(obj$x)] <- -1
###         res <- .C("wholeViterbi", y=as.double(obj$y), x=as.double(obj$x),
###                   genome=as.integer(genome),
###                   index=as.integer(index), k=as.integer(k),
###                   n=as.integer(length(obj$y)), N=as.integer(N), 
###                   mu=as.double(t(obj[[k]]$mu)), sigma.2=as.double(t(obj[[k]]$sigma.2)),
###                   beta=as.double(obj[[k]]$beta),
###                   stat=as.double(obj[[k]]$stat),
###                   filename=as.character(paste(filename, k, sep=""))
###                   )
###       }
###     }
###   }
### }


### ## Example for viterbi
### viterbi.C <- function(y, x=NULL, Chrom=NULL, mu, sigma.2, beta, stat=NULL) {

###   if (!is.null(Chrom)) {
###     index <- diff(Chrom)
###     index <- which(index>0)
###     index <- c(0, index, length(y))
###     genome <- length(index) - 1
###   }
###   else {
###     genome <- 1
###     index <- c(0, length(y))
###   }
###   if (is.null(x)) x <- rep(0, length(y)-1)
###   k <- length(mu)
###   n <- length(y)
###   if (is.null(stat)) stat <- rep(1/k, k)
###   states <- .C("viterbi", y=as.double(y), x=as.double(x),
###                genome=as.integer(genome),
###                index = as.integer(index), k =as.integer(k),
###                n=as.integer(n), mu=as.double(mu),
###                sigma2=as.double(sigma.2), beta=as.double(beta),
###                stat=as.double(stat), states=as.integer(rep(0, n)))
###   states <- states$states
###   states
### }

### ## Maybe we could avoid passing obj
### prob.seq <- function(obj, from, to, filename, alteration="Gain") {

###   if (from > to)
###     stop("'from' must be an integer < than 'to'\n")
###   if (alteration!="Gain" && alteration!="Loss")
###     stop("'alteration' must be either 'Loss' or 'Gain'\n")
###   K <- max(as.numeric(as.character(levels(obj$k))))
###   n <- length(obj$y)
###   probs <- prop.table(table(obj$k))
###   prob.seq <- rep(0, K)
###   for (k in 1:K) {
###     if(probs[k] > 0) {
###       inds <- grep(paste("[", substr(alteration, 1, 1), "]"),
###                    colnames(obj[[k]]$state.labels))
###       if (length(inds) > 0) {
###         seq.iter <- readLines(paste(filename, k, sep=""))
###         seq.iter <- strsplit(seq.iter, "\t")
###         for (i in 1:length(seq.iter)) {
###           seq.iter[[i]] <- as.numeric(seq.iter[[i]])
###           breakpoints <- seq.iter[[i]][seq(from=2, by=2,
###                                            length=length(seq.iter[[i]])/2-1)]
###           breakstates <- seq.iter[[i]][seq(from=1, by=2,
###                                            length=length(seq.iter[[i]])/2-1)]

###           breaks <- c(which.max(from <= breakpoints),
###                       which.max(to <=breakpoints))
### ##           if(is.na(breaks[1]) | is.na(breaks[2])) browser()
###           break.states <- breakstates[breaks[1]:breaks[2]]

###           if (sum(break.states %in% inds)==length(break.states)) {
###             prob.seq[k] <- prob.seq[k] + seq.iter[[i]][length(seq.iter[[i]])]
###           }
###         }

###         prob.seq[k] <- prob.seq[k] / nrow(obj[[k]]$mu)
###       }
###     }
###   }
###   prob.seq <- as.numeric(prob.seq %*% probs)
###   prob.seq
### }  





### getEdges <- function(obj) {
###   K <- max(as.numeric(as.character(levels(obj$k))))
###   probs <- prop.table(table(obj$k))
###   res.tot <- rep(0, length(obj$y))
###   N.tot <- 0
###   for (k in 1:K) {
###       if(probs[k] > 0) {
###           ##       inds <- grep(paste("[", substr(alteration, 1, 1), "]"),
###           ##                    obj[[k]]$state.labels)
###           ##       if (length(inds) > 0) {
###           ## Viterbi of all iterations
###           N <- nrow(obj[[k]]$mu)
###           if (inherits(obj, "RJaCGH.genome")) {
###               index <- diff(obj$Chrom)
###               index <- which(index>0)
###               index <- c(0, index, length(obj$y))
###               genome <- length(index) -1
###           }
###           else {
###               index <- c(0, length(obj$y))
###               genome <- 1
###           }
###           obj$x[is.na(obj$x)] <- -1
###           res <- .C("edges", y=as.double(obj$y), x=as.double(obj$x),
###                     genome=as.integer(genome),
###                     index=as.integer(index), k=as.integer(k),
###                     n=as.integer(length(obj$y)), N=as.integer(N), 
###                     mu=as.double(t(obj[[k]]$mu)),
###                     sigma.2=as.double(t(obj[[k]]$sigma.2)),
###                     beta=as.double(obj[[k]]$beta),
###                     stat=as.double(obj[[k]]$stat),
###                     count_edge = as.integer(rep(0,length(obj$y)))
###                   )
###           res.tot <- res.tot + res$count_edge
###           N.tot <- N.tot + N
###           }
###   }
###   tmp <- res.tot/N.tot
###   if(any(tmp > 1)) stop("a count larger than 1 !!!!")
###   return(tmp)
### }






### readRO1 <- function(x) {
###     fname <- file.path(x$tmpdir, x$fname)
###     if(!file.exists(fname)) {
###         stop(paste("Cannot acces file", fname,
###                    ". Either you are in the wrond directory",
###                    "or the file has been deleted (the later is ",
###                    "a non-recoverable error)."))
###     }
   
###     nameobj <- load(fname)
###     tmp <- get(nameobj)
###     rm(nameobj)
###     return(tmp)
### }










### O <- list()
### O$tmpdir <- "t1"
### O$fname <- "sx3.RData"
### unix.time(cucu <- readRO(O)) ##47.65
### > gc()
###             used  (Mb) gc trigger  (Mb)  max used  (Mb)
### Ncells    134675   7.2     350000  18.7    350000  18.7
### Vcells 100141052 764.1  110743003 845.0 100142729 764.1


### O <- list()
### O$tmpdir <- "t1"
### O$fname <- "sx3.RData"
### unix.time(cucu <- readRO1(O)) ##45.234
### > gc()
###             used  (Mb) gc trigger  (Mb)  max used  (Mb)
### Ncells    134776   7.2     350000  18.7    350000  18.7
### Vcells 100141064 764.1  110743003 845.0 100143380 764.1


### ## time and object size and memory of both read and write

### writeRO <- function(x, rootn, tmpdir) {
###     fn <- tempfile2(rootn, tmpdir, ftype = "RData")
###     .saveRDS(x, file = fn, compress = FALSE)
### }





### ### difference betwen .saveRDS and serialize
### ### zz <- matrix(rnorm(10000 * 1000), ncol = 1000)

### ### rm(ff)
### ### rm(fn)
### ### system("rm sx9.RData")
### ### ff <- "sx9.RData"
### ### unix.time(.saveRDS(zz, ff, compress = FALSE)) # 1, 1, 1

### ### rm(ff)
### ### rm(fn)
### ### system("rm sx9.RData")
### ### ff <- "sx9.RData"
### ### unix.time(.saveRDS(zz, ff)) ## 8, 8, 8


### ### rm(ff)
### ### rm(fn)
### ### system("rm sx9.RData")
### ### ff <- "sx9.RData"
### ### fn <- file(ff, "wb")  ## it has to be binary; ow. toooooo slow
### ### unix.time(serialize(zz, fn)) ## ~ 8
### ### close(fn)








### ### .saveRDS(zz, "t1/sx10.RData", compress = FALSE)
### ### save(file= "t1/sx11.RData", zz, compress = FALSE)

### ### O1 <- list()
### ### O1$tmpdir <- "t1"
### ### O1$fname <- "sx11.RData"
### ### O2 <- list()
### ### O2$tmpdir <- "t1"
### ### O2$fname <- "sx10.RData"

### ### unix.time(cucuA <- readRO2(O1)) ## 1.345
### ### gc()
### ###            used (Mb) gc trigger (Mb) max used (Mb)
### ### Ncells   134848  7.3     350000 18.7   350000 18.7
### ### Vcells 10141116 77.4   11518085 87.9 10142793 77.4


### ### unix.time(cucuB <- readRO(O2)) ## 1.346
### ### gc()
### ### > gc()
### ###            used  (Mb) gc trigger  (Mb) max used  (Mb)
### ### Ncells   134857   7.3     350000  18.7   350000  18.7
### ### Vcells 20141118 153.7   22543147 172.0 20141568 153.7


### ### ## bail out and repeat


### ### unix.time(cucuB <- readRO(O2)) ## 1.346
### ### gc()

### ### unix.time(cucuA <- readRO2(O1)) ## 1.345
### ### gc()

### ### unix.time(save(file= "t1/sx11A.RData", cucuA, compress = FALSE))
### ### unix.time(save(file= "t1/sx11B.RData", cucuA, compress = TRUE))

















## ## iassign.comp <- function(name, obj = obj, res = res) {
## ##     eval(substitute(obj$var <- res$var,
## ##                     list(var = name, obj = obj, res = res)))
## ## }

## assign.components <- function(obj, res) {
##     obj$normal.reference <- normal.reference
##     obj$window <- window
##     obj$singleState <- singleState
##     obj$prob.s <- res$prob.s
##     return(obj)
## }


## assign.components <- function(obj, res) {
##     obj$prob.s <- res$prob.s
##     return(obj)
## }







## ### ### difference betwen .saveRDS and serialize
## nr <- nc <- 3000
## zz <- matrix(rnorm(nr * nc), ncol = nc)

## rm(zzz)
## rm(ff)
## rm(fn)
## system("rm sx9.RData")
## ff <- "sx9.RData"
## unix.time(.saveRDS(zz, ff, compress = FALSE)) # 1, 1, 1
## unix.time(zzz <- .readRDS(ff)) # 1, 1, 1


## rm(ff)
## rm(fn)
## system("rm sx9.RData")
## ff <- "sx9.RData"
## system.time(.saveRDS(zz, ff)) ## 8, 8, 8

## rm(zzz)
## rm(ff)
## rm(fn)
## system("rm sx9.RData")
## ff <- "sx9.RData"
## fn <- file(ff, "wb")
## system.time(serialize(zz, fn)) ##33, 33, 33
## close(fn)

## fn <- file(ff, "rb")
## system.time(zzz <- unserialize(fn))
## close(fn)




## rm(ff)
## rm(fn)
## system("rm sx9.RData")
## ff <- "sx9.RData"
## fn <- file(ff, "w")
## system.time(serialize(zz, fn)) ##33, 33, 33
## close(fn)



