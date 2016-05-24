#******************************************************************************* 
#
# Estimation for Multivariate Normal Data with Monotone Missingness
# Copyright (C) 2007, University of Cambridge
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
#******************************************************************************

##
## bmonomvn.read.traces:
##
## read the traces contained in the files written by the bmonomvn
## C-side, process them as appropriate, and then delete the trace files
##

"bmonomvn.read.traces" <-
  function(N, n, M, nao, oo, nam, capm, mprior, R, cl,
           hs, thin, verb, rmfiles=TRUE)
{
  ## initialize the list
  trace <- list()
  if(verb >= 1) cat("\nGathering traces\n")

  ## read the traces for mu, S
  muS <- read.muS.trace(n, M, nao, oo, nam, verb, rmfiles)
  trace$mu <- muS$mu; trace$S <- muS$S;

  ## read the traces for Data Augmentation
  trace$DA <- read.DA.trace(nao, verb, rmfiles)

  ## read the blasso regression traces
  for(i in 1:length(n)) {
    fname <- paste("blasso_M", i-1, "_n", n[i], ".trace", sep="")
    lname <- paste("M", i-1, ".n", n[i], sep="")
    table <- read.table(fname, header=TRUE)
    trace$reg[[lname]] <- table2blasso(table, thin, mprior, capm,
                                       i-1, n[i], hs, cl)
    if(rmfiles) unlink(fname)

    ## progress meter
    if(verb >= 1) {
      if(i==length(n)) cat(" reg traces 100% done  \r")
      else cat(paste(" reg traces ", round(100*i/length(n)),
                     "% done   \r", sep=""))
    }
  }

  ## cap off with a final newline
  if(verb >= 1) cat("\n")

  return(trace)
}


## read.muS.trace:
##
## read the trace files for mu, S, (.trace)

read.muS.trace <- function(n, M, nao, oo, nam, verb, rmfiles)
  {
    ## initialize the trace list
    trace <- list()
    
    ## read trace of the mean samples (mu)
    if(file.exists(paste("./", "mu.trace", sep=""))) {
      trace$mu <- read.table("mu.trace")
      if(!is.null(oo)) trace$mu <- trace$mu[,oo]
      if(!is.null(nam)) names(trace$mu) <- nam
      if(rmfiles) unlink("mu.trace")
      if(verb >= 1) cat("  mu traces done\n")
    }
    
    ## read trace of the Covar samples (S)
    if(file.exists(paste("./", "S.trace", sep=""))) {
      trace$S <- read.table("S.trace")
      
      ## reorder the columns
      if(!is.null(nao)) {
        om <- matrix(NA, M, M)
        om[lower.tri(om, diag=TRUE)] <- 1:((M+1)*M/2)
        om[upper.tri(om)] <- t(om)[upper.tri(t(om))]
        om <- om[oo,oo]
        om <- om[lower.tri(om, diag=TRUE)]
        trace$S <- trace$S[,om]
      }
      
      ## assign names to the columns
      if(is.null(nam)) nam <- 1:length(n)
      nams <- rep(NA, (M+1)*M/2)
      k <- 1
      for(i in 1:M) for(j in i:M) {
        nams[k] <- paste(nam[i], ":", nam[j], sep="")
        k <- k+1
      }
      names(trace$S) <- nams
      
      ## delete the trace file
      if(rmfiles) unlink("S.trace")
      if(verb >= 1) cat("   S traces done\n")
    }

    ## return the traces
    return(trace)
 }


## read.DA.trace:
##
## read the trace file of the Data Augmentation samples

read.DA.trace <- function(nao, verb, rmfiles)
  {
    DA <- NULL
    
    ## read in the Data Augmentation trace
    if(file.exists(paste("./", "DA.trace", sep=""))) {
      DA <- read.table("DA.trace", header=TRUE)
      
      ## re-name the columns according to nao
      nams <- names(DA)
      for(i in 1:length(nams)) {
        spl <- strsplit(nams[i], "j")
        j <- as.numeric(spl[[1]][2])
        nams[i] <- paste(spl[[1]][1], "j", nao[j], sep="")
      }
      names(DA) <- nams
      
      ## delete the trace file
      if(rmfiles) unlink("DA.trace")
      if(verb >= 1) cat("  DA traces done\n")
    }
    
    ## return traces
    return(DA)
  }



## table2blasso:
##
## change the table trace read in and convert it into a
## skeleton blasso class object so that the blasso methods
## like print, plot, and summary can be used

table2blasso <- function(table, thin, mprior, capm, m, n, hs, cl)
  {
    ## first convert to a list
    tl <- as.list(table)
    
    ## start with the easy scalars
    l <- list(lpost=tl[["lpost"]], llik=tl[["llik"]],
              llik.norm=tl[["llik.norm"]], s2=tl[["s2"]],
              mu=tl[["mu"]], m=tl[["m"]], lambda2=tl[["lambda2"]],
              gamma=tl[["gamma"]], nu=tl[["nu"]], pi=tl[["pi"]])
    
    ## now the vectors
    bi <- grep("beta.[0-9]+", names(table))
    l$beta <- as.matrix(table[,bi])
    ti <- grep("tau2i.[0-9]+", names(table))
    l$tau2i <- as.matrix(table[,ti])
    ## could be that lasso is off, and tau2i doesn't exist
    if(length(l$tau2i) == 0) l$tau2i <- NULL
    else l$tau2i[l$tau2i == -1] <- NA
    ti <- grep("omega2.[0-9]+", names(table))
    l$omega2 <- as.matrix(table[,ti])
    ## could be that Student-t is off, and omega2 doesn't exist
    if(length(l$omega2) == 0) l$omega2 <- NULL

    ## determine horseshoe indicator
    if(!is.null(l$lambda2) && hs) l$hs <- TRUE
    else l$hs <- FALSE

    ## null
    if(is.null(l$llik.norm)) l$llik.norm <- NULL
    if(m == 0) l$beta <- NULL
    
    ## assign "inputs"
    l$T <- length(l$mu)
    l$thin <- "dynamic"
    l$RJ <- !is.null(l$m)
    if(l$RJ) l$mprior <- mprior
    else { l$mprior <- l$m <- l$pi <- NULL }
      
    if(capm) l$M <- min(m, n) 
    else l$M <- m
    
    ## assign the call and the class
    l$call <- cl
    class(l) <- "blasso"
    
    return(l)
  }
