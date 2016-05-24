#' @title A function for computing the spectral information divergence between spectra (sid) 
#' @description
#' This function computes the spectral information divergence (distance) between spectra based on the kullback-leibler divergence algorithm (see details).
#' @usage 
#' sid(Xr, X2 = NULL, 
#'     mode = "density", 
#'     center = FALSE, scaled = TRUE, 
#'     kernel = "gaussian", 
#'     n = if(mode == "density") round(0.5 * ncol(Xr)), 
#'     bw = "nrd0",
#'     reg = 1e-04,
#'     ...)
#' @param Xr a \code{matrix} (or \code{data.frame}) containing the spectral (reference) data.
#' @param X2 an optional \code{matrix} (or \code{data.frame}) containing the spectral data of a second set of samples.
#' @param mode the method to be used for computing the spectral information divergence. Options are \code{"density"} (default) for computing the divergence values on the density distributions of the spectral observations, and \code{"feature"} for computing the divergence vales on the spectral features. See details. 
#' @param center a logical indicating if the computations must be carried out on the centered \code{X} and \code{X2} (if specified) matrices. If \code{mode = "feature"} centering is not carried out since this option does not accept negative values which are generated after centering the matrices. Default is FALSE. See details.
#' @param scaled a logical indicating if the computations must be carried out on the variance scaled \code{X} and \code{X2} (if specified) matrices. Default is TRUE.
#' @param kernel if \code{mode = "density"} a character string indicating the smoothing kernel to be used. It must be one of \code{"gaussian"} (default), \code{"rectangular"}, \code{"triangular"}, \code{"epanechnikov"}, \code{"biweight"}, \code{"cosine"} or \code{"optcosine"}. See the \code{\link[stats]{density}} function of the \code{stats} package.
#' @param n if \code{mode = "density"} a numerical value indicating the number of equally spaced points at which the density is to be estimated. See the \code{\link[stats]{density}} function of the \code{stats} package for further details. Default is \code{round(0.5 * ncol(X))}.
#' @param bw if \code{mode = "density"} a numerical value indicating the smoothing kernel bandwidth to be used. Optionally the character string \code{"nrd0"} can be used, it computes the bandwidth using the \code{\link[stats]{bw.nrd0}} function of the \code{stats} package. See the \code{\link[stats]{density}} and the \code{\link[stats]{bw.nrd0}} functions for more details. By default \code{"nrd0"} is used, in this case the bandwidth is computed as \code{bw.nrd0(as.vector(X))}, if \code{X2} is specified the bandwidth is computed as \code{bw.nrd0(as.vector(rbind(X, X2)))}.
#' @param reg a numerical value higher than 0 which indicates a regularization parameter. Values (probabilities) below this threshold are replaced by this value for numerical stability. Default is 1e-4.
#' @param ... additional arguments to be passed to the \code{\link[stats]{density}} function of the base package. 
#' @details
#' This function computes the spectral information divergence (distance) between spectra. 
#' When \code{mode = "density"}, the function first computes the probability distribution of each spectrum which result in a matrix of density distribution estimates. The density distributions of all the samples in the datasets are compared based on the kullback-leibler divergence algorithm.
#' When \code{mode = "feature"}, the kullback-leibler divergence between all the samples is computed directly on the spectral variables.
#' The spectral information divergence (SID) algorithm (Chang, 2000) uses the Kullback-Leibler divergence (\eqn{KL}) or relative entropy (Kullback and Leibler, 1951) to account for the vis-NIR information provided by each spectrum. The SID between two spectra (\eqn{x_{i}} and \eqn{x_{j}}) is computed as follows:
#' \deqn{
#'    SID(x_{i},x_{j}) = KL(x_{i} \left |\right | x_{j}) + KL(x_{j} \left |\right | x_{i})
#'    }  
#' \deqn{
#'    SID(x_{i},x_{j}) = \sum_{l=1}^{k} p_l \ log(\frac{p_l}{q_l}) + \sum_{l=1}^{k} q_l \ log(\frac{q_l}{p_l})
#' }
#' where \eqn{k} represents the number of variables or spectral features, \eqn{p} and \eqn{q} are the probability vectors of \eqn{x_{i}} and \eqn{x_{j}} respectively which are calculated as:
#' \deqn{
#'    p = \frac{x_i}{\sum_{l=1}^{k} x_{i,l}}
#'    }
#' \deqn{   
#'    q = \frac{x_j}{\sum_{l=1}^{k} x_{j,l}}
#'    }
#' From the above equations it can be seen that the original SID algorithm assumes that all the components in the data matrices are nonnegative. Therefore centering cannot be applied when \code{mode = "feature"}. If a data matrix with negative values is provided and \code{mode = "feature"}, the \code{sid} function automatically scales the matrix as follows:
#' \deqn{
#'     X_s = \frac{X-min(X)}{max(X)-min(X)}
#'    }
#' or 
#' \deqn{
#' X_{s} = \frac{X-min(X, X2)}{max(X, X2)-min(X, X2)}
#'    }
#' \deqn{
#' X2_{s} = \frac{X2-min(X, X2)}{max(X, X2)-min(X, X2)}
#'    }
#' if \code{X2} is specified. The 0 values are replaced by a regularization parameter (\code{reg} argument) for numerical stability.
#' The default of the \code{sid} function is to compute the SID based on the density distributions of the spectra (\code{mode = "density"}). For each spectrum in \code{X} the density distribution is computed using the \code{\link[stats]{density}} function of the \code{stats} package.
#' The 0 values of the estimated density distributions of the spectra are replaced by a regularization parameter (\code{"reg"} argument) for numerical stability. Finally the divergence between the computed spectral histogramas is computed using the SID algorithm. Note that if \code{mode = "density"}, the \code{sid} function will accept negative values and matrix centering will be possible.
#' @return a \code{list} with the following components:
#' \itemize{
#'  \item{\code{sid}}{ if only \code{"X"} is specified (i.e. \code{X2 = NULL}), a square symmetric matrix of SID distances between all the components in \code{"X"}. If both \code{"X"} and \code{"X2"} are specified, a matrix of SID distances between the components in \code{"X"} and the components in \code{"X2"}) where the rows represent the objects in \code{"X"} and the columns represent the objects in \code{"X2"}}
#'  \item{\code{Xr}}{ the (centered and/or scaled if specified) spectral \code{X} matrix}
#'  \item{\code{X2}}{ the (centered and/or scaled if specified) spectral \code{X2} matrix}
#'  \item{\code{densityDisXr}}{ if \code{mode = "density"}, the computed density distributions of \code{Xr}}
#'  \item{\code{densityDisX2}}{ if \code{mode = "density"}, the computed density distributions of \code{X2}}
#'  }
#' @references Chang, C.I. 2000. An information theoretic-based approach to spectral variability, similarity and discriminability for hyperspectral image analysis. IEEE Transactions on Information Theory 46, 1927-1932.
#' @seealso \code{\link[stats]{density}}
#' @author Leonardo Ramirez-Lopez
#' @examples
#' \dontrun{
#' require(prospectr)
#' 
#' data(NIRsoil)
#' 
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' Xu <- Xu[!is.na(Yu),]
#' Xr <- Xr[!is.na(Yr),]
#' 
#' # Example 1
#' # Compute the SID distance between all the samples in Xr
#' xr.sid <- sid(Xr = Xr)
#' xr.sid
#' 
#' # Example 2
#' # Compute the SID distance between the samples in Xr and the samples 
#' # in Xu
#' xru.sid <- sid(Xr = Xr, X2 = Xu)
#' xru.sid
#' 
#' # Example 3
#' # Compute the SID distance between the samples in Xr and the samples 
#' # in Xu using the histograms
#' xru.sid.hist <- sid(Xr = Xr, X2 = Xu, mode = "feature")
#' xru.sid.hist
#' }
#' @export

#######################################################################
# resemble
# Copyrigth (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#######################################################################

## History:
## 09.03.2014 Leo     Doc examples  were formated with a max. line width
##                       


sid <- function(Xr, X2 = NULL, 
                mode = "density", 
                center = FALSE, scaled = TRUE, 
                kernel = "gaussian", 
                n = if(mode == "density") round(0.5 * ncol(Xr)), 
                bw = "nrd0", 
                reg = 1e-4, ...){
  
  # Sanity checks
  if(!is.logical(center))
    stop("'center' argument must be logical")
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical")
  
  if(!is.numeric(reg))
    stop("'reg' argument must be numeric")
  
  if(length(reg) > 1)
    stop("'reg' must be a single numeric value")
  
  n <- n
  X <- Xr
  rm(Xr)
  if(is.null(X2)){
    if(mode == "density"){
      
      match.arg(kernel, c("gaussian", "rectangular", "triangular", "biweight", "epanechnikov", "cosine", "optcosine"))
      
      if(center){
        X <- sweep(x = X, MARGIN = 2, FUN = "-", STATS = colMeans(X))
      }
      if(scaled){      
        X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = apply(X, 2, sd))
      }
      Xo <- X
      
      if(bw == "nrd0")
        bw <- bw.nrd0(as.vector(X))
      
      rng <- c(0.65 * min(X), 1.35 * max(X))
      #rng <- c(0.80 * min(X), 1.2 * max(X))
      str <- seq(rng[1], rng[2], diff(rng)/(n-1))
      
      densM <- matrix(0, nrow(X), length(str))
      colnames(densM) <- str
      for(i in 1:nrow(X))
      {
        densM[i,] <- density(X[i,], from = rng[1], to = rng[2], n = n, kernel = kernel, bw = bw, ...)$y
      }
      X <- densM
      
      reg <- 10^-4
      w <- X < reg
      X[w] <- reg
      rm(densM)
      
      ps <- sweep(x = (X), MARGIN = 1, FUN = "/", STATS = rowSums((X)))
      qs <- sweep(x = (X), MARGIN = 1, FUN = "/", STATS = rowSums((X)))
      
      
      abKld <- matrix(NA, nrow(X), nrow(X))
      for(i in 1:nrow(X))
      {
        dvsn <- sweep(x = ps, MARGIN = 2, FUN = "/", STATS = qs[i,])
        abKld[i,] <- rowSums(ps * (log(dvsn)))
      }
      klKl <- (abKld + t(abKld))/2
      colnames(klKl) <- paste("Xr", 1:nrow(X), sep = ".")
      rownames(klKl) <- paste("Xr", 1:nrow(X), sep = ".")
      rownames(X) <- rownames(klKl) 
      return(list(sid = klKl, Xr = Xo, densityDisXr = X, xdval = str))
    }
    
    if(mode == "feature"){
      if(center){
        warning("centering was not applied since when mode = 'feature', the function does not accept negative values (which are generated by centering)")
      }
      if(sum(abs(X)) != (sum(X))){
        X <- (X - min(X))/diff(range(X))
        message(paste("the 'X' matrix was rescaled from 0 to 1 values, since it contained negative values.", if(scaled){"After this the variables were scaled to unit variance."}))
      }
      if(scaled){      
        X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = apply(X, 2, sd))
      }
      Xo <- X
      
      reg <- 10^-4
      w <- X < reg
      X[w] <- reg
      
      ps <- sweep(x = (X), MARGIN = 1, FUN = "/", STATS = rowSums((X)))
      qs <- sweep(x = (X), MARGIN = 1, FUN = "/", STATS = rowSums((X)))
      
      
      abKld <- matrix(NA, nrow(X), nrow(X))
      for(i in 1:nrow(X))
      {
        dvsn <- sweep(x = ps, MARGIN = 2, FUN = "/", STATS = qs[i,])
        abKld[i,] <- rowSums(ps * (log(dvsn)))
      }
      klKl <- (abKld + t(abKld))/2
      colnames(klKl) <- paste("Xr", 1:nrow(X), sep = ".")
      rownames(klKl) <- paste("Xr", 1:nrow(X), sep = ".")
      return(list(sid = klKl, Xr = X))
    }
  }
  
  
  if(!is.null(X2)){
    if(ncol(X) != ncol(X2))
      stop("The number of variables in X must be equal to the number of variables in X2")
    
    if(mode == "density"){
      
      match.arg(kernel, c("gaussian", "rectangular", "triangular", "biweight", "epanechnikov", "cosine", "optcosine"))
      
      X <- rbind(X, X2)
      if(center){
        X <- sweep(x = X, MARGIN = 2, FUN = "-", STATS = colMeans(X))
      }
      if(scaled){      
        X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = apply(X, 2, sd))
      }
      Xo <- X
      
      
      if(bw == "nrd0")
        bw <- bw.nrd0(as.vector(X))
      
      rng <- c(0.65 * min(X), 1.35 * max(X))
      #rng <- c(0.80 * min(X), 1.2 * max(X))
      str <- seq(rng[1], rng[2], diff(rng)/(n-1))
      
      densM <- matrix(0, nrow(X), length(str))
      colnames(densM) <- str
      for(i in 1:nrow(X))
      {
        densM[i,] <- density(X[i,], from = rng[1], to = rng[2], n = n, kernel = kernel, bw = bw, ...)$y
      }
      
      reg <- 10^-4
      w <- densM < reg
      densM[w] <- reg
      
      X2 <- densM[(nrow(X) - nrow(X2) +1):nrow(densM), ,drop = FALSE]
      X <- densM[1:(nrow(densM) - nrow(X2)), ,drop = FALSE]
      rm(densM)
      
      ps <- sweep(x = X, MARGIN = 1, FUN="/", STATS = rowSums(X))
      qs <- sweep(x = X2, MARGIN = 1, FUN="/", STATS = rowSums(X2))
      
      klKl <- matrix(NA, nrow(X), nrow(X2))
      for(i in 1:nrow(X2))
      {
        dvsnA <- sweep(x = ps, MARGIN = 2, FUN = "/", STATS = qs[i,])
        dvsnB <- sweep(x = ps^-1, MARGIN = 2, FUN = "/", STATS = (qs[i,]^-1))
        klKl[,i] <-  (rowSums(ps * log(dvsnA)) + colSums(qs[i,] * t(log(dvsnB))))/2
      }
      colnames(klKl) <- paste("X2", 1:nrow(X2), sep = ".")
      rownames(klKl) <- paste("Xr", 1:nrow(X), sep = ".")
      rownames(X) <- rownames(klKl)
      rownames(X2) <- colnames(klKl) 
      
      return(list(sid = klKl, Xr = Xo[1:nrow(X), ,drop = FALSE], X2 = Xo[(nrow(Xo) - nrow(X2) + 1):nrow(Xo), ,drop = FALSE], 
                  densityDisXr = X, densityDisX2 = X2, 
                  xdval = str))               
    }
    
    if(mode == "feature"){
      X <- rbind(X, X2)
      if(center){
        warning("centering was not applied since when mode = 'feature', the function does not accept negative values (which are generated by centering)")
      }
      if(sum(abs(X)) != (sum(X))){
        X <- (X - min(X))/diff(range(X))
        message(paste("since either or both 'X' and 'X2' contained negative values, these matrices were rescaled  to a 0 to 1 range (see details on the 'sid' function).", if(scaled){"After this the variables were scaled to unit variance."}, "Finally the matrix was split back into 'X' and 'X2'."))
      }
      if(scaled){      
        X <- sweep(x = X, MARGIN = 2, FUN = "/", STATS = apply(X, 2, sd))
      }
      
      X2 <- X[(nrow(X) - nrow(X2) + 1):nrow(X), ]
      X <- X[1:(nrow(X) - nrow(X2)), ]
      
      ps <- sweep(x = X, MARGIN = 1, FUN="/", STATS = rowSums(X))
      qs <- sweep(x = X2, MARGIN = 1, FUN="/", STATS = rowSums(X2))
      
      klKl <- matrix(NA, nrow(X), nrow(X2))
      for(i in 1:nrow(X2))
      {
        dvsnA <- sweep(x = ps, MARGIN = 2, FUN = "/", STATS = qs[i,])
        dvsnB <- sweep(x = ps^-1, MARGIN = 2, FUN = "/", STATS = (qs[i,]^-1))
        klKl[,i] <-  (rowSums(ps * log(dvsnA)) + colSums(qs[i,] * t(log(dvsnB))))/2
      }
      
      colnames(klKl) <- paste("X2", 1:nrow(X2), sep = ".")
      rownames(klKl) <- paste("Xr", 1:nrow(X), sep = ".")
      
      return(list(sid = klKl, Xr = X, X2 = X2))
    }
  }
}
