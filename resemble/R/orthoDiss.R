#' @title A function for computing dissimilarity matrices from orthogonal projections (orthoDiss)
#' @description
#' This function computes dissimilarities (in an orthogonal space) between either observations in a given set or between observations in two different sets.
#' The dissimilarities are computed based on either principal component projection or partial least squares projection of the data. After projecting the data, 
#' the Mahalanobis distance is applied.
#' @usage 
#' orthoDiss(Xr, X2 = NULL, 
#'           Yr = NULL, 
#'           pcSelection = list("cumvar", 0.99), 
#'           method = "pca", 
#'           local = FALSE, 
#'           k0, 
#'           center = TRUE, scaled = FALSE, 
#'           return.all = FALSE, cores = 1, ...)
#' @param Xr a \code{matrix} (or \code{data.frame}) containing the (reference) data.
#' @param X2 an optional \code{matrix} (or \code{data.frame}) containing data of a second set of observations(samples).
#' @param Yr either if the method used in the \code{pcSelection} argument is \code{"opc"} or if the \code{sm} argument is either \code{"pls"} or \code{"loc.pls"}, then it must be a \code{vector} containing the side information corresponding to the spectra in \code{Xr}. It is equivalent to the \code{sideInf} parameter of the \code{\link{simEval}} function. It can be a numeric \code{vector} or \code{matrix} (regarding one or more continuous variables). The root mean square of differences (rmsd) is used for assessing the similarity between the samples and their corresponding most similar samples in terms of the side information provided. When \code{sm = "pc"}, this parameter can also be a single discrete variable of class \code{factor}. In such a case the kappa index is used. See \code{\link{simEval}} function for more details.
#' @param pcSelection a list which specifies the method to be used for identifying the number of principal components to be retained for computing the Mahalanobis distance of each sample in \code{sm = "Xu"} to the centre of \code{sm = "Xr"}. It also specifies the number of components in any of the following cases: \code{sm = "pc"}, \code{sm = "loc.pc"}, \code{sm = "pls"} and \code{sm = "loc.pls"}. This list must contain two objects in the following order: \itemize{
#'        \item{\code{method}:}{the method for selecting the number of components. Possible options are:  \code{"opc"} (optimized pc selection based on Ramirez-Lopez et al. (2013a, 2013b). See the \code{\link{orthoProjection}} function for more details;  \code{"cumvar"} (for selecting the number of principal components based on a given cumulative amount of explained variance); \code{"var"} (for selecting the number of principal components based on a given amount of explained variance); and  \code{"manual"} (for specifying manually the desired number of principal components)}
#'        \item{\code{value}:}{a numerical value that complements the selected method. If \code{"opc"} is chosen, it must be a value indicating the maximal number of principal components to be tested (see Ramirez-Lopez et al., 2013a, 2013b). If \code{"cumvar"} is chosen, it must be a value (higher than 0 and lower than 1) indicating the maximum amount of cumulative variance that the retained components should explain. If \code{"var"} is chosen, it must be a value (higher than 0 and lower than 1) indicating that components that explain (individually) a variance lower than this threshold must be excluded. If \code{"manual"} is chosen, it must be a value specifying the desired number of principal components to retain.
#'        }}
#'        The default method for the \code{pcSelection} argument is \code{"opc"} and the maximal number of principal components to be tested is set to 40.
#'        Optionally, the \code{pcSelection} argument admits \code{"opc"} or \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character string. In such a case the default for \code{"value"} when either \code{"opc"} or \code{"manual"} are used is 40. When \code{"cumvar"} is used the default \code{"value"} is set to 0.99 and when \code{"var"} is used the default \code{"value"} is set to 0.01.
#' @param method the method for projecting the data. Options are: "pca" (principal component analysis using the singular value decomposition algorithm), "pca.nipals" (principal component analysis using the non-linear iterative partial least squares algorithm) and "pls" (partial least squares). See the \code{\link{orthoProjection}} function for further details on the projection methods.
#' @param local a logical indicating whether or not to compute the distances locally (i.e. projecting locally the data) by using the \eqn{k0} nearest neighbour samples of each sample. Default is \code{FALSE}. See details.
#' @param k0 if \code{local = TRUE} a numeric integer value which indicates the number of nearest neighbours(\eqn{k0}) to retain in order to recompute the local orthogonal distances.
#' @param center a logical indicating if the spectral data \code{Xr} (and \code{X2} if specified) must be centered. If \code{X2} is specified the data is centered on the basis of \eqn{Xr \cup Xu}. For dissimilarity computations based on pls, the data is always centered for the projections. 
#' @param scaled a logical indicating if \code{Xr} (and \code{X2} if specified) must be scaled. If \code{X2} is specified the data is scaled on the basis of \eqn{Xr \cup Xu}.
#' @param return.all a logical. In case \code{X2} is specified it indicates whether or not the distances between all the elements resulting from \eqn{Xr \cup Xu} must be computed.
#' @param cores number of cores used when \code{method} in \code{pcSelection} is \code{"opc"} (which can be computationally intensive) and \code{local = FALSE} (default = 1). Dee details.
#' @param ... additional arguments to be passed to the \code{\link{orthoProjection}} function.
#' @details
#' When \code{local = TRUE}, first a global distance matrix is computed based on the parameters specified. Then, by using this matrix for each target observation, a given set of nearest neighbours (\eqn{k0}) are identified. These neighbours (together with the target observation) are projected (from the original data space) onto a (local) orthogonal space (using the same parameters specified in the function). 
#' In this projected space the Mahalanobis distance between the target sample and the neighbours is recomputed. A missing value is assigned to the samples that do not belong to this set of neighbours (non-neighbour samples).
#' In this case the dissimilarity matrix cannot be considered as a distance metric since it does not necessarily satisfies the symmetry condition for distance matrices (i.e. given two samples \eqn{x_i} and \eqn{x_j}, the local dissimilarity (\eqn{d}) between them is relative since generally \eqn{d(x_i, x_j) \neq d(x_j, x_i)}). On the other hand, when \code{local = FALSE}, the dissimilarity matrix obtained can be considered as a distance matrix.
#' @return a \code{list} of class \code{orthoDiss} with the following components:
#' \itemize{
#'  \item{\code{n.components}}{ the number of components (either principal components or partial least squares components) used for computing the global distances.}
#'  \item{\code{global.variance.info}}{ the information about the expalined variance(s) of the projection. When \code{local = TRUE}, the information corresponds to the global projection done prior computing the local projections.}
#'  \item{\code{loc.n.components}}{ if \code{local = TRUE}, a \code{data.frame} which specifies the number of local components (either principal components or partial least squares components) used for computing the dissimilarity between each target sample and its neighbour samples.}
#'  \item{\code{dissimilarity}}{ the computed dissimilarity matrix. If \code{local = FALSE} a distance \code{matrix}. If \code{local = TRUE} a \code{matrix} of class \code{orthoDiss}. In this case each column represent the dissimilarity between a target sample and its neighbourhood.}
#'  }
#' Multi-threading for the computation of dissimilarities (see \code{cores} parameter) is based on OpenMP and hence works only on windows and linux.
#' @author Leonardo Ramirez-Lopez
#' @references 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local approach for modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196, 268-279.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{orthoProjection}}, \code{\link{simEval}}
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
#' Yu <- Yu[!is.na(Yu)]
#' 
#' Xr <- Xr[!is.na(Yr),]
#' Yr <- Yr[!is.na(Yr)] 
#' 
#' # Computation of the orthogonal dissimilarity matrix using the 
#' # default parameters
#' ex1 <- orthoDiss(Xr = Xr, X2 = Xu)
#' 
#' # Computation of a principal component dissimilarity matrix using 
#' # the "opc" method for the selection of the principal components
#' ex2 <- orthoDiss(Xr = Xr, X2 = Xu, 
#'                  Yr = Yr, 
#'                  pcSelection = list("opc", 40), 
#'                  method = "pca", 
#'                  return.all = TRUE)
#' 
#' # Computation of a partial least squares (PLS) dissimilarity 
#' # matrix using the "opc" method for the selection of the PLS 
#' # components
#' ex3 <- orthoDiss(Xr = Xr, X2 = Xu, 
#'                  Yr = Yr, 
#'                  pcSelection = list("opc", 40), 
#'                  method = "pls")
#' 
#' # Computation of a partial least squares (PLS) local dissimilarity 
#' # matrix using the "opc" method for the selection of the PLS 
#' # components
#' ex4 <- orthoDiss(Xr = Xr, X2 = Xu, 
#'                  Yr = Yr, 
#'                  pcSelection = list("opc", 40), 
#'                  method = "pls",
#'                  local = TRUE,
#'                  k0 = 200)
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
## 09.03.2014 Leo     In the doc was specified that multi-threading is 
##                    not working for mac
## 13.03.2014 Antoine The explanation of the cores argument was modified
## 07.09.2014 Antoine A bug handling Yr as a matrix was fixed
## 02.12.2015 Leo     The function now outputs an the object global.variance.info
##                    which provides information on the explained variance
##                    of the projections.

orthoDiss <- function(Xr, X2 = NULL, 
                      Yr = NULL, 
                      pcSelection = list("cumvar", 0.99), 
                      method = "pca", 
                      local = FALSE, 
                      k0, 
                      center = TRUE, scaled = FALSE, 
                      return.all = FALSE, cores = 1, ...){
   
  in.call <- match.call()
  if(!is.null(in.call$call.))
    call. <- in.call$call.
  else
    call. <- TRUE
  
  if(!is.numeric(cores))
    stop("The 'cores' argument must be numeric")
  
  if(!is.logical(local))
    stop("'local' argument must be logical", call. = call.)
  
  if(!is.logical(return.all))
    stop("'return.all' argument must be logical", call. = call.)
  
  if(!is.logical(center))
    stop("'center' argument must be logical", call. = call.)
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical", call. = call.)
  
  if(local){
    if(k0 < 10)
      stop("k0 cannot be smaller than 10", call. = call.)
    if(k0 > nrow(Xr))
      stop("k0 cannot be larger than the numer of elements in Xr", call. = call.)
  }
  
  if(!is.null(Yr))
  {
    if(is.null(dim(Yr)))
      Yr <- matrix(Yr, length(Yr))
    Yr <- as.data.frame(Yr)
  }else{
    if(pcSelection[[1]] == "opc" | method == "pls"){
      stop("Yu must be provided either when the 'opc' is used in pcSelection is used or method = 'pls'", call. = call.)
    }
  }
  
  prj <- orthoProjection(Xr = Xr, Yr = Yr, X2 = X2, method = method, pcSelection = pcSelection, center = center, scaled = scaled, cores = cores, call. = FALSE)
  scores <- prj$scores
  scores <- sweep(prj$scores, 2, prj$sc.sdv, "/")
  n.components <- prj$n.components
  if(is.null(X2)){
    distnc <- fDiss(Xr = scores, X2 = NULL, method = "euclid", center = FALSE, scaled = FALSE)
    dimnames(distnc) <- list(rownames(scores), rownames(scores))
  }else{
    if(!return.all)
    {
      distnc <- fDiss(Xr = scores[1:nrow(Xr), ,drop = FALSE], X2 = scores[(1+nrow(Xr)):nrow(scores), ,drop = FALSE], method = "euclid", center = FALSE, scaled = FALSE)
      dimnames(distnc) <- list(rownames(scores[1:nrow(Xr),]), rownames(scores[(1+nrow(Xr)):nrow(scores),]))
    }else{
      distnc <- fDiss(Xr = scores, X2 = NULL, method = "euclid", center = FALSE, scaled = FALSE)
      dimnames(distnc) <- list(rownames(scores), rownames(scores))
    }
  }
  
  d.val <- i <- NULL
  
  if(local)
  {
    
    if(is.null(X2)){
      
      srch <- foreach(d.val = iter(distnc, by = "column"), i = icount(),
                      .combine = cbind, .inorder = FALSE,
                      .export = c("orthoProjection", "pcProjection", "plsProjection", "simEval"), .packages = c("resemble")) %dopar% { # Parallel computations
                        
                        sel <- order(d.val)[1:k0]
                        
                        pca.i <- orthoProjection(Xr = Xr[sel,], Yr = Yr[sel,], X2 = NULL, method = method, pcSelection = pcSelection, center = center, scaled = scaled, call. = FALSE)
                        scores <- sweep(pca.i$scores, 2, pca.i$sc.sdv, "/")
                        
                        dst.i <- data.frame(c(NaN, sel), c(pca.i$n.components, 0, fDiss(Xr = scores[1,,drop = FALSE], X2 = scores[-1,,drop = FALSE], method = "euclid", center = FALSE, scaled = FALSE)))
                        names(dst.i) <-  c("sel", i)
                        dst.i
                      }
      
      ind <- rep(as.numeric(srch[2, seq(1,ncol(srch),2)]), each = 2) 
      ind[seq(1,length(ind),2)] <- ind[seq(1,length(ind),2)] + 0.1
      ind[seq(2,length(ind),2)] <- ind[seq(2,length(ind),2)] + 0.2
      ind <- order(ind)
      srch <- srch[,ind]
      loc.n.components <- as.numeric(srch[1, seq(2,length(ind),2)])
      srch <- srch[-1,]
      nms <- dimnames(distnc) 
      distnc <- matrix(NaN, nrow(Xr), nrow(Xr))
      dimnames(distnc) <- nms
      count <- 0
      for(j in seq(1,length(ind),2))
      {
        count <- count + 1
        distnc[srch[,j],count] <- srch[,j+1]
      }
      colnames(distnc) <- paste("Xr", 1:ncol(distnc), sep = ".")
      rownames(distnc) <- colnames(distnc)
    }else{
      if(!return.all)
      {
        srch <- foreach(d.val = iter(distnc, by = "column"), i = icount(),
                        .combine = cbind, .inorder = FALSE,
                        .export = c("orthoProjection", "pcProjection", "plsProjection", "simEval"), .packages = c("resemble")) %dopar% { # Parallel computations
                          
                          sel <- order(d.val)[1:k0]
                          
                          pca.i <- orthoProjection(Xr = Xr[sel,], Yr = Yr[sel,,drop = FALSE], X2 = X2[i, ,drop = FALSE], method = method, pcSelection = pcSelection, center = center, scaled = scaled, call. = FALSE)
                          scores <- sweep(pca.i$scores, 2, pca.i$sc.sdv, "/")
                          dst.i <- data.frame(c(NaN, i, sel), c(pca.i$n.components, 0, fDiss(Xr = scores[nrow(scores), ,drop=FALSE], X2 = scores[-nrow(scores),,drop=FALSE], method = "euclid", center = FALSE, scaled = FALSE)))
                          names(dst.i) <-  c("sel", i)
                          dst.i
                        }
        
        ind <- rep(as.numeric(srch[2, seq(1,ncol(srch),2)]), each = 2) 
        ind[seq(1,length(ind),2)] <- ind[seq(1,length(ind),2)] + 0.1
        ind[seq(2,length(ind),2)] <- ind[seq(2,length(ind),2)] + 0.2
        ind <- order(ind)
        srch <- srch[,ind]
        loc.n.components <- as.numeric(srch[1, seq(2,length(ind),2)])
        srch <- srch[-1,]
        nms <- dimnames(distnc) 
        distnc <- matrix(NaN, nrow(Xr), nrow(X2))
        dimnames(distnc) <- nms
        count <- 0
        for(j in seq(1,length(ind),2))
        {
          count <- count + 1
          distnc[srch[-1,j],count] <- srch[-1,j+1]
        }
        colnames(distnc) <- paste("X2", 1:ncol(distnc), sep = ".")
        rownames(distnc) <- paste("Xr", 1:nrow(distnc), sep = ".")
      }else{
        if("opc" %in% pcSelection | method == "pls")
        {
          rowOrd <- (apply(distnc[,(1 + nrow(Xr)):ncol(distnc)], 2, order))[1:k0,]
          cnt <- colSums(rowOrd <= nrow(Xr))
          miss <- paste((1:nrow(X2))[cnt < 10])
          
          if(length(miss) > 0)
          {
            if(length(miss) > 80)
              miss <- c(miss[1:80], "...")
            warning(ngettext(length(miss), "variable", "The samples in X2 with the following indexes can be problematic:\n"),
                    paste(sQuote(miss), collapse = ", "),
                    ngettext(length(miss), "contains", " \n Most of their k0 neighbours are in X2 and less than 10 neighbours are in Xr \n therefore the side information (Yu) can be insuficient to produce reliable \n estimations of the optimal number of principal components"))
          }
        }
        nr <- nrow(Xr)
        nu <- nrow(X2)
        
        X <- rbind(Xr, X2)
        Y <- rbind(Yr, rep(NaN, nrow(X2)))
        rm(Xr)
        rm(X2)
        
        srch <- foreach(d.val = iter(distnc, by = "column"), i = icount(),
                        .combine = cbind, .inorder = FALSE,
                        .export = c("simEval"), .packages = c("resemble")) %dopar% { # Parallel computations
                          
                          sel <- order(d.val)[1:k0]
                          pca.i <- prcomp(X[sel,], center = center, scale = scaled)
                          
                          pca.i <- orthoProjection(Xr = X[sel,], Yr = Y[sel,,drop = FALSE], X2 = NULL, method = method, pcSelection = pcSelection, center = center, scaled = scaled, call. = FALSE)
                          scores <- sweep(pca.i$scores, 2, pca.i$sc.sdv, "/")
                          dst.i <- data.frame(c(NaN, sel), c(pca.i$n.components, 0, fDiss(Xr = scores[1, ,drop=FALSE], X2 = scores[-1,,drop=FALSE], method = "euclid", center = FALSE, scaled = FALSE)))
                          names(dst.i) <-  c("sel", i)
                          dst.i
                        }
        
        ind <- rep(as.numeric(srch[2, seq(1,ncol(srch),2)]), each = 2) 
        ind[seq(1,length(ind),2)] <- ind[seq(1,length(ind),2)] + 0.1
        ind[seq(2,length(ind),2)] <- ind[seq(2,length(ind),2)] + 0.2
        ind <- order(ind)
        srch <- srch[,ind]
        loc.n.components <- as.numeric(srch[1, seq(2,length(ind),2)])
        srch <- srch[-1,]
        nms <- dimnames(distnc) 
        distnc <- matrix(NaN, nrow(X), nrow(X))
        dimnames(distnc) <- nms
        count <- 0
        for(j in seq(1,length(ind),2))
        {
          count <- count + 1
          distnc[srch[,j],count] <- srch[,j+1]
        }
        colnames(distnc) <- c(paste("Xr", 1:nr, sep = "."), paste("X2", 1:nu, sep = "."))
        rownames(distnc) <- colnames(distnc)
        
      }
    }
    resultsList <- list(n.components = n.components, global.variance.info = prj$variance, loc.n.components = data.frame(sample.nm = colnames(distnc), sample = 1:ncol(distnc), loc.n.components = loc.n.components), dissimilarity = distnc)
    class(resultsList) <- c("orthoDiss", "list") 
    class(resultsList$dissimilarity) <- c("localOrthoDiss","matrix")
    return(resultsList)    
  }else{
    resultsList <- list(n.components = n.components, global.variance.info = prj$variance, dissimilarity = distnc)
    class(resultsList) <- c("orthoDiss", "list") 
    class(resultsList$dissimilarity) <- c("orthoDiss","matrix") 
    return(resultsList)
  }
}
