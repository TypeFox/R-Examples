#' @title Orthogonal projections using partial least squares and principal component analysis
#' @aliases orthoProjection 
#' @aliases plsProjection 
#' @aliases pcProjection 
#' @aliases predict.orthoProjection
#' @description
#' Functions to perform orthogonal projections of high dimensional data matrices using partial least squares (pls) and principal component analysis (pca)
#' @usage 
#' orthoProjection(Xr, X2 = NULL, 
#'                 Yr = NULL, 
#'                 method = "pca", pcSelection = list("cumvar", 0.99), 
#'                 center = TRUE, scaled = FALSE, cores = 1, ...)
#'                 
#' pcProjection(Xr, X2 = NULL, Yr = NULL, 
#'              pcSelection = list("cumvar", 0.99), 
#'              center = TRUE, scaled = FALSE, 
#'              method = "pca",
#'              tol = 1e-6, max.iter = 1000, 
#'              cores = 1, ...)  
#'               
#' plsProjection(Xr, X2 = NULL, Yr, 
#'               pcSelection = list("opc", 40), 
#'               scaled = FALSE, 
#'               tol = 1e-6, max.iter = 1000, 
#'               cores = 1, ...) 
#'               
#' \method{predict}{orthoProjection}(object, newdata, ...)
#'                            
#' @param Xr a \code{matrix} (or \code{data.frame}) containing the (reference) data.
#' @param X2 an optional \code{matrix} (or \code{data.frame}) containing data of a second set of observations(samples).
#' @param Yr if the method used in the \code{pcSelection} argument is \code{"opc"} or if the \code{sm} argument is either \code{"pls"} or \code{"loc.pls"}, then it must be a \code{vector} containing the side information corresponding to the spectra in \code{Xr}. It is equivalent to the \code{sideInf} parameter of the \code{\link{simEval}} function. It can be a numeric \code{vector} or \code{matrix} (regarding one or more continuous variables). The root mean square of differences (rmsd) is used for assessing the similarity between the samples and their corresponding most similar samples in terms of the side information provided. When \code{sm = "pc"}, this parameter can also be a single discrete variable of class \code{factor}. In such a case the kappa index is used. See \code{\link{simEval}} function for more details.
#' @param pcSelection a list which specifies the method to be used for identifying the number of principal components to be retained for computing the Mahalanobis distance of each sample in \code{sm = "Xu"} to the centre of \code{sm = "Xr"}. It also specifies the number of components in any of the following cases: \code{sm = "pc"}, \code{sm = "loc.pc"}, \code{sm = "pls"} and \code{sm = "loc.pls"}. This list must contain two objects in the following order: \itemize{
#'        \item{\code{method}:}{the method for selecting the number of components. Possible options are:  \code{"opc"} (optimized pc selection based on Ramirez-Lopez et al. (2013a, 2013b) in which the side information concept is used, see details), \code{"cumvar"} (for selecting the number of principal components based on a given cumulative amount of explained variance); \code{"var"} (for selecting the number of principal components based on a given amount of explained variance); and  \code{"manual"} (for specifying manually the desired number of principal components)}
#'        \item{\code{value}:}{a numerical value that complements the selected method. If \code{"opc"} is chosen, it must be a value indicating the maximal number of principal components to be tested (see Ramirez-Lopez et al., 2013a, 2013b). If \code{"cumvar"} is chosen, it must be a value (higher than 0 and lower than 1) indicating the maximum amount of cumulative variance that the retained components should explain. If \code{"var"} is chosen, it must be a value (higher than 0 and lower than 1) indicating that components that explain (individually) a variance lower than this threshold must be excluded. If \code{"manual"} is chosen, it must be a value specifying the desired number of principal components to retain.
#'        }}
#'        The default method for the \code{pcSelection} argument is \code{"opc"} and the maximal number of principal components to be tested is set to 40.
#'        Optionally, the \code{pcSelection} argument admits \code{"opc"} or \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character string. In such a case the default for \code{"value"} when either \code{"opc"} or \code{"manual"} are used is 40. When \code{"cumvar"} is used the default \code{"value"} is set to 0.99 and when \code{"var"} is used the default \code{"value"} is set to 0.01.
#' @param method the method for projecting the data. Options are: "pca" (principal component analysis using the singular value decomposition algorithm), "pca.nipals" (principal component analysis using the non-linear iterative partial least squares algorithm) and "pls" (partial least squares).
#' @param center a logical indicating if the data \code{Xr} (and \code{X2} if specified) must be centered. If \code{X2} is specified the data is centered on the basis of \eqn{Xr \cup Xu}. This argument only applies to the principal components projection. For pls projections the data is always centered. 
#' @param scaled a logical indicating if \code{Xr} (and \code{X2} if specified) must be scaled. If \code{X2} is specified the data is scaled on the basis of \eqn{Xr \cup Xu}.
#' @param tol tolerance limit for convergence of the algorithm in the nipals algorithm (default is 1e-06). In the case of PLS this applies only to Yr with more than two variables.
#' @param max.iter maximum number of iterations (default is 1000). In the case of \code{method = "pls"} this applies only to \code{Yr} matrices with more than one variable.
#' @param cores number of cores used when \code{method} in \code{pcSelection} is \code{"opc"} (which can be computationally intensive) (default = 1). Dee details.
#' @param ... additional arguments to be passed to \code{pcProjection} or \code{plsProjection}.
#' @param object  object of class "orthoProjection" (as returned by \code{orthoProjection}, \code{pcProjection} or \code{plsProjection}).
#' @param newdata an optional data frame or matrix in which to look for variables with which to predict. If omitted, the scores are used. It must contain the same number of columns, to be used in the same order.
#' @details
#' In the case of \code{method = "pca"}, the algrithm used is the singular value decomposition in which given a data matrix \eqn{X}, is factorized as follows:
#' \deqn{
#'      X = UDV^{\mathrm{T}}
#'      }
#' where \eqn{U} and \eqn{V} are othogonal matrices, and where \eqn{U} is a matrix of the left singular vectors of \eqn{X}, \eqn{D} is a diagonal matrix containing the singular values of \eqn{X} and \eqn{V} is the is a matrix of the right singular vectors of \eqn{X}.
#' The matrix of principal component scores is obtained by a matrix multiplication of \eqn{U} and \eqn{D}, and the matrix of principal component loadings is equivalent to the matrix \eqn{V}. 
#' When \code{method = "pca.nipals"}, the algorithm used for principal component analysis is the non-linear iterative partial least squares (nipals).
#' In the case of the of the partial least squares projection (a.k.a projection to latent structures) the nipals regression algorithm. Details on the "nipals" algorithm are presented in Martens (1991).
#' When \code{method = "opc"}, the selection of the components is carried out by using an iterative method based on the side information concept (Ramirez-Lopez et al. 2013a, 2013b). First let be \eqn{P} a sequence of retained components (so that \eqn{P = 1, 2, ...,k }. 
#' At each iteration, the function computes a dissimilarity matrix retaining \eqn{p_i} components. The values of the side information of the samples are compared against the side information values of their most spectrally similar samples. 
#' The optimal number of components retrieved by the function is the one that minimizes the root mean squared differences (RMSD) in the case of continuous variables, or maximizes the kappa index in the case of categorical variables. In this process the \code{\link{simEval}} function is used. 
#' Note that for the \code{"opc"} method is necessary to specify \code{Yr} (the side information of the samples).
#' Multi-threading for the computation of dissimilarities (see \code{cores} parameter) is based on OpenMP and hence works only on windows and linux. 
#' @return \code{orthoProjection}, \code{pcProjection}, \code{plsProjection}, return a \code{list} of class \code{orthoProjection} with the following components:
#' \itemize{
#'  \item{\code{scores}}{ a \code{matrix} of scores corresponding to the samples in \code{Xr} and \code{X2} (if it applies). The number of components that the scores represent is given by the number of components chosen in the function.}
#'  \item{\code{X.loadings}}{ a \code{matrix} of loadings corresponding to the explanatory variables. The number of components that these loadings represent is given by the number of components chosen in the function.}
#'  \item{\code{Y.loadings}}{ a \code{matrix} of partial least squares loadings corresponding to \code{Yr}. The number of components that these loadings represent is given by the number of components chosen in the function. This object is only returned if the partial least squares algorithm was used.}
#'  \item{\code{weigths}}{ a \code{matrix} of partial least squares ("pls") weights. This object is only returned if the "pls" algorithm was used.}
#'  \item{\code{projectionM}}{ a \code{matrix} that can be used to project new data onto a "pls" space. This object is only returned if the "pls" algorithm was used.}
#'  \item{\code{variance}}{ a \code{matrix} indicating the standard deviation of each component (sdv), the cumulative explained variance (cumExplVar) and the variance explained by each single component (explVar). These values are computed based on the data used to create the projection matrices. 
#'                            For example if the "pls" method was used, then these values are computed based only on the data that contains information on \code{Yr} (i.e. the \code{Xr} data)
#'                            If the principal component method is used, the this data is computed on the basis of \code{Xr} and \code{X2} (if it applies) since both matrices are employed in the computation of the projection matrix (loadings in this case)}. 
#'  \item{\code{svd}}{ the standard deviation of the retrieved scores.}
#'  \item{\code{n.components}}{ the number of components (either principal components or partial least squares components) used for computing the global distances.}
#'  \item{\code{opcEval}}{ a \code{data.frame} containing the statistics computed for optimizing the number of principal components based on the variable(s) specified in the \code{Yr} argument. If \code{Yr} was a continuous  was a continuous \code{vector} or \code{matrix} then this object indicates the root mean square of differences (rmse) for each number of components. If \code{Yr} was a categorical variable this object indicates the kappa values for each number of components. 
#'                        This object is returned only if \code{"opc"} was used within the \code{pcSelection} argument. See the \code{\link{simEval}} function for more details.}
#'  \item{\code{method}}{ the \code{orthoProjection} method used.}
#'  }
#'  \code{predict.orthoProjection}, returns a matrix of scores proprojected for \code{newdtata}.
#' @author Leonardo Ramirez-Lopez
#' @references 
#' Martens, H. (1991). Multivariate calibration. John Wiley & Sons.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local approach for modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196, 268-279.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{orthoDiss}}, \code{\link{simEval}}, \code{\link{mbl}}
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
#' # A partial least squares projection using the "opc" method
#' # for the selection of the optimal number of components
#' plsProj <- orthoProjection(Xr = Xr, Yr = Yr, X2 = Xu, 
#'                            method = "pls", 
#'                            pcSelection = list("opc", 40))
#'                            
#' # A principal components projection using the "opc" method
#' # for the selection of the optimal number of components
#' pcProj <- orthoProjection(Xr = Xr, Yr = Yr, X2 = Xu, 
#'                           method = "pca", 
#'                           pcSelection = list("opc", 40))
#'                            
#' # A partial least squares projection using the "cumvar" method
#' # for the selection of the optimal number of components
#' plsProj2 <- orthoProjection(Xr = Xr, Yr = Yr, X2 = Xu, 
#'                             method = "pls", 
#'                             pcSelection = list("cumvar", 0.99))
#' } 
#' @rdname orthoProjection                           
#' @export

######################################################################
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
######################################################################

## History:
## 09.03.2014 Leo     In the doc was specified that multi-threading is 
##                    not working for mac
## 10.03.2014 Leo     Sanity check for the number of cases in Yr and Xr
## 13.03.2014 Antoine The explanation of the cores argument was modified
## 03.12.2015 Leo     The pls projection function is now based on Rcpp
## 04.12.2015 Leo     The center argument was removed from the pls projection 
##                    function. Matrices are now always centered.
## 04.12.2015 Leo     the dist function needs to be replaced by the fDiss. dist() retrieves
##                    squared distances!
## 07.12.2015 Leo     The dist function was removed from these functions (it was implemented)
##                    for testing for a while but it was not finally released.
## 07.12.2015 Leo     The pcProjection, plsProjection and the predict.othoProjection are now
##                    visible.

orthoProjection <- function(Xr, X2 = NULL, 
                            Yr = NULL, 
                            method = "pca", pcSelection = list("cumvar", 0.99), 
                            center = TRUE, scaled = FALSE, 
                            cores = 1, ...)
{
  in.call <- match.call()
  if(!is.null(in.call$call.))
    call. <- in.call$call.
  else
    call. <- TRUE
  
  method <- match.arg(method, c("pls", "pca", "pca.nipals"))
  
  if(method == "pls")
  {
    if(!is.numeric(as.matrix(Yr)))
      stop("When pls projection is used, 'Yr' must be numeric", call. = call.)
    proj <- plsProjection(Xr = Xr, Yr = Yr, X2 = X2, pcSelection = pcSelection, scaled = scaled, cores = cores, call. = FALSE, ...)    
    mthd <- "pls (nipals)"
  }else{
    mthd <- ifelse(method == "pca", "pca (svd)", "pca (nipals)")
    proj <- pcProjection(Xr = Xr, Yr = Yr, X2 = X2, method , pcSelection = pcSelection, 
                         center = center, scaled = scaled, 
                         method = method, cores = cores, call. = FALSE, ...)
  }
  gc() # free some memory 
  proj$method <- mthd
  class(proj) <- c("orthoProjection","list")  
  return(proj)                      
}


#' @rdname orthoProjection      
#' @aliases orthoProjection 
#' @aliases plsProjection 
#' @aliases pcProjection 
#' @aliases predict.orthoProjection
#' @export
pcProjection <- function(Xr, X2 = NULL, Yr = NULL, 
                         pcSelection = list("cumvar", 0.99), 
                         center = TRUE, scaled = FALSE, 
                         method = "pca",
                         tol = 1e-6, max.iter = 1000, 
                         cores = 1, ...){
  in.call <- match.call()
  if(!is.null(in.call$call.))
    call. <- in.call$call.
  else
    call. <- TRUE
  
  pcSel <- pcSelection[[1]]
  match.arg(pcSel, c("opc", "var", "cumvar", "manual"))
  
  match.arg(method, c("pca", "pca.nipals"))
  
  if(!is.logical(center))
    stop("'center' must be logical", call. = call.)
  
  if(!is.logical(scaled))
    stop("'scaled' must be logical", call. = call.)
  
  if(is.factor(Yr))
    if(length(Yr) != nrow(Xr))
      stop("The number of obsevations in 'Yr' does not match with the number of observations in 'Xr'. In addition, when 'Yr' is a factor, only one variable is accepted", call. = call.)
  
  nxr <- nrow(Xr)
  
  if(!is.null(X2))
  {
    if(pcSel %in% c("opc", "manual"))
    {
      if(length(pcSelection) == 1){
        trsh <- min(nrow(Xr) + nrow(X2), ncol(Xr))
        trsh <- ifelse(trsh > 40, 40, trsh)
        pcSelection <- list(method = pcSel, value = trsh)
        message(paste("Since the value of the 'pcSelection' argument is missing, it was automatically set to ", trsh, sep=""))
      } else {
        if(!is.list(pcSelection))
          stop("The pcSelection argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifiying only the method can be used, in this case the parameter value is set automatically", call. = call.)
        pcSelection <- list(method = pcSelection[[1]], value = floor(pcSelection[[2]]))
        if(!is.numeric(pcSelection$value)) 
          stop("The second object in 'pcSelection' must be an integer value", call. = call.)
        if(pcSelection$value < 2 | pcSelection$value > min(nrow(Xr) + nrow(X2), ncol(Xr)))
          stop(paste("The maximum number of principal components must be a value between 2 and", min(nrow(Xr) + nrow(X2), ncol(Xr))), call. = call.)
      }
      max.i <- pcSelection$value
    }
    nr <- nrow(Xr)
    n2 <- nrow(X2)
    Xr <- rbind(Xr, X2)
  }else{
    if(pcSel %in% c("opc", "manual"))
    {
      if(length(pcSelection) == 1){
        trsh <- min(dim(Xr))
        trsh <- ifelse(trsh > 40, 40, trsh)
        pcSelection <- list(method = pcSel, value = trsh)
        message(paste("Since the value of the 'pcSelection' argument is missing, it was automatically set to ", trsh, sep=""))
      } else {
        if(!is.list(pcSelection))
          stop("The pcSelection argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifiying only the method can be used, in this case the parameter value is set automatically", call. = call.)     
        pcSelection <- list(method = pcSelection[[1]], value = floor(pcSelection[[2]]))
        if(!is.numeric(pcSelection$value)) 
          stop("The second object in 'pcSelection' must be an integer value", call. = call.)
        if(pcSelection$value < 2 | pcSelection$value > min(dim(Xr)))
          stop(paste("The maximum number of principal components must be a value between 2 and", min(dim(Xr))), call. = call.)
      }
      max.i <- pcSelection$value
    }
    nr <- nrow(Xr)
    n2 <- NULL
  }
  
  if(pcSel %in% c("cumvar", "var"))
  {
    if(length(pcSelection) == 1){
      if(pcSel == "cumvar"){
        pcSelection <- list("cumvar", 0.99)
        message(paste("Since the value of the 'pcSelection' argument is missing, the amount of cumulative variance that the components to be retained should explain was automatically set to 0.99 (99%)"))
      }else{
        pcSelection <- list("var", 0.01)
        message(paste("Since the value of the 'pcSelection' argument is missing, the amount of variance that the last component to be retained should explain was automatically set to 0.01 (1%)"))
      }
      names(pcSelection) <- c("method", "value")
    } else {
      if(!is.list(pcSelection))
        stop("The 'pcSelection' argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifiying only the method can be used, in this case the parameter value is set automatically", call. = call.)     
      pcSelection <- list(pcSelection[[1]], pcSelection[[2]])
      names(pcSelection) <- c("method", "value")
      if(!is.numeric(pcSelection$value)) 
        stop("The second object in 'pcSelection' must be a numeric value", call. = call.)
      if(pcSelection$value > 1 | pcSelection$value <= 0) 
        stop(paste("When the method for 'pcSelection' is either 'var' or 'cumvar' the value in 'pcSelection' must be a number higher than 0 and lower than/or equal to 1"), call. = call.)
    }
    max.i <- min(dim(Xr)) - 1
  }
  
  psel <- pcSelection
  
  # center 
  if(center){
    cvec <- colMeans(Xr)
    X0 <- sweep(x = Xr, MARGIN = 2, FUN = "-", STATS = cvec)
  }else{
    cvec <- rep(0, ncol(Xr))
    X0 <- Xr
  }
  
  if(scaled)
  {
    sf <- cSds(X0)
    X0 <- sweep(x = X0, MARGIN = 2, FUN = "/", STATS = sf)
  } else{
    sf <- rep(1, ncol(X0))
  }
  
  if(method == "pca")
  {
    svDecomp <- svd(x = X0)
    # Loadings and scores
    pcLoadings <- t(svDecomp$v)
    pcScores <- svDecomp$u %*% diag(svDecomp$d)
    # assign names
    colnames(pcScores) <- paste("pc", 1:ncol(pcScores), sep = "")
    rownames(pcScores) <- c(paste("Xr.", 1:nr, sep = ""), if(!is.null(X2)){paste("X2.", 1:n2, sep = "")})
    colnames(pcLoadings) <- colnames(X0)
    rownames(pcLoadings) <- paste("pc", 1:nrow(pcLoadings), sep = "")
    
    # Variance of each PC variable
    sdPC <- (cSds(pcScores))
    # Compute the percentage of explained variance for all the PCs
    ons <- (svDecomp$d)^2 / (nrow(X0) - 1)
    ev.i <- ons / sum(ons)
    ev <- diffinv(ev.i[-1], xi = ev.i[1])
    variance <- rbind(sdv = sdPC,  cumExplVar = ev, explVar = ev.i)  
    
    if(pcSel == "var"){
      ev <- 1 - ev.i
      pcSelection$value <- 1 - pcSelection$value
    }
  }
  
  
  if(method == "pca.nipals")
  {
    xvar <- sum(cSds(X0)^2)
    nPf <- min(dim(X0))
    pcScores <- matrix(NA, nrow(X0), nPf)
    pcLoadings <- matrix(NA, nPf, ncol(X0))
    exv <- matrix(NA, 3, nPf)
    
    colnames(pcScores) <- paste("pc", 1:ncol(pcScores), sep = "")
    rownames(pcScores) <- c(paste("Xr.", 1:nr, sep = ""), if(!is.null(X2)){paste("X2.", 1:n2, sep = "")})
    colnames(pcLoadings) <- colnames(X0)
    rownames(pcLoadings) <- paste("pc", 1:nrow(pcLoadings), sep = "")
    
    if(pcSel %in% c("opc", "manual"))
      pcSelection$value <- pcSelection$value - 1
    
    
    xx <- X0
    for(i in 1:max.i)
    {
      tt <- xx[,1]
      iter <- keepg <- T
      while(keepg){
        pp <- (t(xx) %*% tt)/(t(tt) %*% tt)[[1]]
        pp <- pp / ((t(pp)%*%pp)^0.5)[[1]]
        #tt.i <- (t(pp)%*%t(xx) / (t(pp)%*%pp)[[1]]) # ths can be removed (t(pp)%*%pp)[[1]]
        tt.i <- t(pp)%*%t(xx)
        val.n <- (((tt-tt.i)%*%t(tt-tt.i)))^0.5
        keepg <- !(val.n <= tol)
        if(max.iter <= iter)
          keepg <- FALSE
        iter <- 1 + iter
        tt <- tt.i[1,]
      }
      xx <- xx - (tt %*% t(pp))
      pcScores[,i]<- tt
      xpvar <- sum(cSds(xx)^2)
      pcLoadings[i,] <- pp[,1]
      exv[1,i] <- (sd(pcScores[,i]))^2
      ev <- (sum(exv[1,1:i]) / sum(xvar))
      exv[2,i] <- ev 
      exv[3,i] <- exv[1,i]/sum(xvar)
      
      if(pcSel %in% c("var", "cumvar"))
      {
        if(ifelse(pcSel == "cumvar", ev > pcSelection$value, exv[3,i] < pcSelection$value))
        {
          i <- i - 1
          if(i == 0)
            stop("With the current value in the 'pcSelection' argument, no components are selected. Try another value.", call. = call.)
          break
        }
      }  
    }
    sel <- i
    pcScores <- pcScores[ ,1:sel , drop=FALSE]
    pcLoadings <- pcLoadings[1:sel, , drop=FALSE]
    ev.i <- exv[3,]
    variance <- rbind(sdv = exv[1,],  cumExplVar = exv[3,], explVar = exv[2,]) 
    colnames(variance) <- paste("pc", 1:ncol(variance), sep = "")
  }
  
  
  if(pcSel == "opc")
  {
    if(missing(Yr) | is.null(Yr))
      stop("'Yr' must be provided when the 'opc' method is used for selecting the optimal number of principal components", call. = call.)
    
    if(!is.factor(Yr))
    {
      Yr <- as.matrix(Yr)
      ny <- ncol(Yr)
      if(nrow(Yr) != nxr)
        stop("The number of rows in Xr does not match the number of cases in Yr")
    }else{ 
      ny <- 1
    }
    
    if(ny > 1)
    {
      if(sum(duplicated(colnames(Yr))) > 0)
        stop("names of the Yr variables must be different")
      result <- matrix(NA, pcSelection$value, ny + 1)
      s.scores <- sweep(pcScores, MARGIN = 2, STATS = cSds(pcScores), FUN = "/")      
      
      
      #       for(i in 1:pcSelection$value){        
      #         sc <- s.scores[1:nr, i, drop = FALSE]       
      #         if(i==1){ # this allows to build-up the distance matrix progressively, saving time
      #           if(cores>1)
      #             d <- fastDistVV(sc, cores)
      #           else
      #             d <- dist(sc)^2 # The square is necessary here to get the results compatible with fastDistVV
      #         } else {
      #           if(cores>1)
      #             d <- d + fastDistVV(sc, cores)
      #           else
      #             d <- d + dist(sc)^2 # The square is necessary here to get the results compatible with fastDistVV
      #         }
      #         tmp <- simEval(d = d, sideInf = Yr, call. = FALSE, lower.tri = TRUE,  cores = cores)        
      #         result[i,1:ny] <- getElement(tmp$eval, "rmsd")
      #         result[i,1 + ny] <- getElement(tmp$global.eval, "mn.sd.rmsd")
      #       }
      
      # this builds-up the distance matrices progressively, saving time
      
      if(.Platform$OS.type == "unix" & cores > 1) {
        for (i in 1:pcSelection$value) {
          sc <- s.scores[1:nr, i, drop = FALSE]
          if (i == 1) {
            d <- fastDistVV(sc, cores)
          }
          else {
            d <- d + fastDistVV(sc, cores)
          }
          tmp <- simEval(d = d, sideInf = Yr, call. = FALSE, 
                         lower.tri = TRUE, cores = cores)
          result[i, 1:ny] <- getElement(tmp$eval, "rmsd")
          result[i, 1 + ny] <- getElement(tmp$global.eval, 
                                          "mn.sd.rmsd")
        }
      }else{
        for (i in 1:pcSelection$value) {
          sc <- s.scores[1:nr, i, drop = FALSE]
          if (i == 1) {
            d <- fastDistVV(sc, cores)
          }
          else {
            d <- d + fastDistVVL(sc)
          }
          tmp <- simEval(d = d, sideInf = Yr, call. = FALSE, 
                         lower.tri = TRUE, cores = cores)
          result[i, 1:ny] <- getElement(tmp$eval, "rmsd")
          result[i, 1 + ny] <- getElement(tmp$global.eval, 
                                          "mn.sd.rmsd")
        }
      }
      
      results <- data.frame(1:pcSelection$value, result)
      colnames(results) <- c("pc", paste("rmsd.Y", 1:ny, sep = ""), "mn.sd.rmsd.Y")
      nPf <- which.min(results$mn.sd.rmsd.Y)
    }else{
      if(!is.factor(Yr))
        ext <- "rmsd"
      else
        ext <- "kappa"      
      
      result <- rep(NA, pcSelection$value)
      s.scores <- sweep(pcScores, MARGIN = 2, STATS = cSds(pcScores), FUN = "/")
      
      
#       for(i in 1:pcSelection$value){           
#         sc <- s.scores[1:nr,i, drop=FALSE]
#         if(i==1){ # this allows to build-up the distance matrix progressively, saving time
#           if(cores>1)
#             d <- fastDistVV(sc, cores)
#           else
#             d <- dist(sc)
#         } else {
#           if(cores>1)
#             d <- d + fastDistVV(sc, cores)
#           else
#             d <- d + dist(sc)
#         }
#         tmp <- simEval(d = d, sideInf = Yr, call. = FALSE, lower.tri = TRUE, cores = cores)
#         result[i] <- getElement(tmp$eval, ext)
#       }
      
      for (i in 1:pcSelection$value) {
        sc <- s.scores[1:nr, i, drop = FALSE]
        if (i == 1) {
          d <- fastDistVV(sc, cores)
        }
        else {
          d <- d + fastDistVV(sc, cores)
        }
        tmp <- simEval(d = d, sideInf = Yr, call. = FALSE, 
                       lower.tri = TRUE, cores = cores)
        result[i] <- getElement(tmp$eval, ext)
      }
      
      
      #result <- result[order(result[,1]),]
      results <- data.frame(1:pcSelection$value, result)
      colnames(results) <- c("pc", names(tmp$eval)[[1]])
      nPf <- which.min(results[,names(tmp$eval)[[1]]])
    }
    sel <- nPf
    sel[1:nPf] <- TRUE
  }
  
  if(pcSel %in% c("cumvar","var") & method == "pca"){
    sel <- (ev <= pcSelection$value)
    sel <- sum(sel)
  }
  
  if(pcSel == "manual" & method == "pca") {
    sel <- (1:ncol(pcScores)) <= pcSelection$value
    sel <- sum(sel)
  }
  
  if(pcSel == "opc") {
    sel <- sum(sel)
  }
  
  fresults <- list(scores = pcScores[,1:sel, drop = FALSE], 
                   X.loadings = pcLoadings[1:sel,, drop = FALSE], 
                   variance = variance[,1:sel, drop = FALSE], 
                   sc.sdv = cSds(pcScores[ ,1:sel, drop = FALSE]), 
                   n.components = sel, pcSelection = psel, 
                   center = cvec,
                   scale = sf)
  colnames(fresults$variance) <- rownames(fresults$X.loadings)
  fresults$method <- ifelse(method == "pca", "pca (svd)", "pca (nipals)")
  if(pcSel == "opc") {fresults$opcEval <- results} 
  class(fresults) <- c("orthoProjection","list") 
  return(fresults)
}

#' @rdname orthoProjection      
#' @aliases orthoProjection 
#' @aliases plsProjection 
#' @aliases pcProjection 
#' @aliases predict.orthoProjection
#' @export
plsProjection <- function(Xr, X2 = NULL, Yr, 
                          pcSelection = list("opc", 40), 
                          scaled = FALSE, 
                          tol = 1e-6, max.iter = 1000, 
                          cores = 1, ...){
  in.call <- match.call()
  if(!is.null(in.call$call.))
    call. <- in.call$call.
  else
    call. <- TRUE
  
  pcSel <- pcSelection[[1]]
  match.arg(pcSel, c("opc", "var", "cumvar", "manual"))
  
  if(!is.numeric(cores))
    stop("The 'cores' argument must be numeric")
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical", call. = call.)
  
  if(missing(Yr))
    stop("'Yr' must be provided", call. = call.)
  
  if(!is.numeric(as.matrix(Yr)))
    stop("'Yr' must be numeric", call. = call.)
  
  if(!is.null(X2))
  {
    if(pcSel %in% c("opc", "manual"))
    {
      if(length(pcSelection) == 1){
        trsh <- min(dim(Xr))
        trsh <- ifelse(trsh > 40, 40, trsh)
        pcSelection <- list(method = pcSel, value = trsh)
        message(paste("Since the value of the 'pcSelection' argument is missing, it was automatically set to ", trsh, sep=""))
      } else {
        if(!is.list(pcSelection))
          stop("The pcSelection argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifiying only the method can be used, in this case the parameter value is set automatically", call. = call.)
        pcSelection <- list(method = pcSelection[[1]], value = floor(pcSelection[[2]]))
        if(!is.numeric(pcSelection$value)) 
          stop("The second object in 'pcSelection' must be an integer value", call. = call.)
        if(pcSelection$value < 2 | pcSelection$value > min(dim(Xr))) 
          stop(paste("The maximum number of components must be a value between 2 and", min(dim(Xr))), call. = call.)
      }
      max.i <- pcSelection$value
    }
  }else{
    if(pcSel %in% c("opc", "manual"))
    {
      if(length(pcSelection) == 1){
        trsh <- min(dim(Xr))
        trsh <- ifelse(trsh > 40, 40, trsh)
        pcSelection <- list(method = pcSel, value = trsh)
        message(paste("Since the value of the 'pcSelection' argument is missing, it was automatically set to ", trsh, sep=""))
      } else {
        if(!is.list(pcSelection))
          stop("The pcSelection argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifiying only the method can be used, in this case the parameter value is set automatically", call. = call.)     
        pcSelection <- list(method = pcSelection[[1]], value = floor(pcSelection[[2]]))
        if(!is.numeric(pcSelection$value)) 
          stop("The second object in 'pcSelection' must be an integer value", call. = call.)
        if(pcSelection$value < 2 | pcSelection$value > min(dim(Xr)))
          stop(paste("The maximum number of components must be a value between 2 and", min(dim(Xr))), call. = call.)  
      }
      max.i <- pcSelection$value
    }
  }
  
  if(pcSel %in% c("cumvar", "var"))
  {
    if(length(pcSelection) == 1){
      if(pcSel == "cumvar"){
        pcSelection <- list("cumvar", 0.99)
        message(paste("Since the value of the 'pcSelection' argument is missing, the amount of cumulative variance that the components to be retained should explain was automatically set to 0.99 (99%)"))
      }else{
        pcSelection <- list("var", 0.01)
        message(paste("Since the value of the 'pcSelection' argument is missing, the amount of variance that the last component to be retained should explain was automatically set to 0.01 (1%)"))
      }
      names(pcSelection) <- c("method", "value")
    } else {
      if(!is.list(pcSelection))
        stop("The 'pcSelection' argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifiying only the method can be used, in this case the parameter value is set automatically", call. = call.)     
      pcSelection <- list(pcSelection[[1]], pcSelection[[2]])
      names(pcSelection) <- c("method", "value")
      if(!is.numeric(pcSelection$value)) 
        stop("The second object in 'pcSelection' must be a numeric value", call. = call.)
      if(pcSelection$value > 1 | pcSelection$value <= 0) 
        stop(paste("When the method for 'pcSelection' is either 'var' or 'cumvar' the value in 'pcSelection' must be a number higher than 0 and lower than/or equal to 1"), call. = call.)
    }
    max.i <- min(dim(Xr)) - 1
  }
  
  psel <- pcSelection
  Yr <- as.matrix(Yr)
  if(nrow(Yr) != nrow(Xr))
    stop("The number of rows in Xr does not match the number of cases in Yr")
  nas <- rowSums(is.na(Yr)) > 0
  
  ny <- ncol(Yr)
  
  X0 <- Xr
  Y0 <- Yr  
  inx.in <- 1:nrow(Xr)
  if(sum(nas) > 0)
  {
    inx.out <- (1:nrow(Xr))[nas]
    inx.in <- (1:nrow(Xr))[!nas]
    Xout <- Xr[inx.out, ]
    X0 <- Xr[inx.in, ]
    Y0 <- Yr[inx.in, , drop = FALSE]
    if(pcSel %in% c("opc", "manual"))
      if(min(dim(Xr)) < pcSelection$value)
        stop("Since there are some missing values in Yr, the corresponding observations (including the Xr observations) cannot be used for pls projection. The problem is that number of components specified in the 'pcSelection' argument exceeds the number of remaining observations. Try another number of components.", call. = call.)
  }
  
  nPf <- min(dim(X0))
  weights <- matrix(NA, nPf, ncol(X0))
  scores <- matrix(NA, nrow(X0), nPf)
  X.loadings <- matrix(NA, nPf, ncol(X0))
  Y.loadings <- matrix(NA, nPf, ny)
  exv <- matrix(NA, 3, nPf)
  
  if(pcSel %in% c("opc", "manual")){
    pcSelection$value <- pcSelection$value - 1
    plsp <- opls(X = X0, 
                 Y = as.matrix(Y0), 
                 ncomp = max.i,
                 scale = scaled,            
                 maxiter = max.iter,
                 tol = tol,
                 regression = TRUE)
    nPf <- plsp$ncomp
  }else{
    plsp <- opls(X = X0, 
                 Y = as.matrix(Y0), 
                 ncomp = max.i,
                 scale = scaled,            
                 maxiter = max.iter,
                 tol = tol,
                 regression = FALSE,
                 pcSelmethod = pcSel,
                 pcSelvalue = pcSelection$value)
    nPf <- plsp$ncomp
  }
  
  ev <- 0
  
  
  
  # If the number of PLS factors is optimized by using the opc method then...
  if(pcSel == "opc")
  {
    pcSelection$value <- pcSelection$value + 1 
    if(ny > 1)
    {
      result <- matrix(NA, pcSelection$value, ny + 1)
      s.scores <- sweep(plsp$scores, MARGIN =2, STATS = cSds(plsp$scores), FUN = "/")
      
      
#       for(i in 1:pcSelection$value){
#         sc <- s.scores[,i, drop=FALSE]
#         if(i==1){ # this allows to build-up the distance matrix progressively, saving time
#           if(cores>1)
#             d <- fastDistVV(sc, cores)
#           else
#             d <- dist(sc)^2 # The square is necessary here to get the results compatible with fastDistVV
#         } else {
#           if(cores>1)
#             d <- d + fastDistVV(sc, cores)
#           else
#             d <- d + (dist(sc)^2) # The square is necessary here to get the results compatible with fastDistVV
#         }
#         tmp <- simEval(d = d, sideInf = Yr[inx.in,,drop=FALSE], lower.tri = TRUE, cores = cores)
#         result[i,1:ny] <- getElement(tmp$eval, "rmsd")
#         result[i,1+ny] <- getElement(tmp$global.eval, "mn.sd.rmsd")
#       }
      
      for (i in 1:pcSelection$value) {
        sc <- s.scores[, i, drop = FALSE]
        if (i == 1) {
          d <- fastDistVV(sc, cores)
        }
        else {
          d <- d + fastDistVV(sc, cores)
        }
        tmp <- simEval(d = d, sideInf = Yr[inx.in, ,drop = FALSE], lower.tri = TRUE, cores = cores)
        result[i, 1:ny] <- getElement(tmp$eval, "rmsd")
        result[i, 1 + ny] <- getElement(tmp$global.eval, "mn.sd.rmsd")
      }
      
      
      results <- data.frame(1:pcSelection$value, result)
      colnames(results) <- c("pls", paste("rmsd.Y", 1:ny, sep = ""), "mn.sd.rmsd.Y")
      nPf <- which.min(results$mn.sd.rmsd.Y)
    }else{
      result <- rep(NA, pcSelection$value)
      s.scores <- sweep(plsp$scores, MARGIN = 2, STATS = cSds(plsp$scores), FUN = "/")
      
      
#       for(i in 1:pcSelection$value){
#         sc <- s.scores[,i, drop = FALSE]
#         if(i == 1){ # this allows to build-up the distance matrix progressively, saving time
#           if(cores>1){
#             d <- fastDistVV(sc, cores)
#           }else
#             d <- dist(sc)^0.5
#         } else {
#           if(cores>1)
#             d <- d + fastDistVV(sc, cores)
#           else
#             d <- d + (dist(sc)^0.5)
#         }
#         tmp <- simEval(d = d, sideInf = Yr[inx.in,,drop=FALSE], lower.tri = TRUE, cores = cores)
#         result[i] <- getElement(tmp$eval, "rmsd")
#       }
#       
      
      
      for (i in 1:pcSelection$value) {
        sc <- s.scores[, i, drop = FALSE]
        if (i == 1) {
          d <- fastDistVV(sc, cores)
        }
        else {
          d <- d + fastDistVV(sc, cores)
        }
        tmp <- simEval(d = d, sideInf = Yr[inx.in, , 
                                           drop = FALSE], lower.tri = TRUE, cores = cores)
        result[i] <- getElement(tmp$eval, "rmsd")
      }
      
      
      
      results <- data.frame(1:pcSelection$value, result)
      colnames(results) <- c("pls", "rmsd.Y")
      nPf <- which.min(results$rmsd.Y)
    }
  }
  
  
  # Select the necessary components
  exv <- plsp$variance$x.var[,1:nPf, drop=FALSE]
  weights <- plsp$weights[1:nPf, , drop=FALSE]
  scores <- plsp$scores[,1:nPf, drop=FALSE]
  X.loadings <- plsp$X.loadings[1:nPf, , drop=FALSE]
  Y.loadings <- plsp$Y.loadings[1:nPf, , drop=FALSE]
  
  # Give some names...
  colnames(X.loadings) <- colnames(X0)
  rownames(X.loadings) <- paste("pls", 1:nPf, sep="")
  rownames(Y.loadings) <- rownames(X.loadings)
  colnames(Y.loadings) <- paste("Y.pls", 1:ny, sep ="")
  rownames(weights) <- rownames(X.loadings)
  colnames(weights) <- colnames(X.loadings)
  colnames(exv) <- rownames(X.loadings)
  rownames(exv) <- c("sdv", "cumExplVar.X", "explVar.X")
  
  yex <- plsp$variance$y.var[,1:nPf,drop = FALSE]
  
  colnames(yex) <- rownames(Y.loadings)
  rownames(yex) <- paste("explVar.Yr", 1:ny, sep = "")
  
  if(sum(nas) > 0)
  {
    scores.a <- matrix(NA, length(c(inx.in, inx.out)), ncol(scores)) 
    scores.a[inx.out,] <- projectpls(projectionm = plsp$projectionM,
                                     ncomp = nPf, 
                                     newdata = Xout,
                                     scale = scaled,
                                     Xcenter = plsp$transf$Xcenter,
                                     Xscale = plsp$transf$Xscale)
    scores.a[inx.in,] <- scores
    scores <- scores.a
  }
  
  rownames(scores) <- paste("Xr.", 1:nrow(scores), sep = "")
  colnames(scores) <- rownames(Y.loadings)
  
  
  if(!is.null(X2)){
    if(is.vector(X2)){
      X2 <- t(X2)
    }
    scoresX2 <- projectpls(projectionm = plsp$projectionM,
                           ncomp = nPf, 
                           newdata = X2,
                           scale = scaled,
                           Xcenter = plsp$transf$Xcenter,
                           Xscale = plsp$transf$Xscale)
    
    colnames(scoresX2) <- rownames(Y.loadings)
    rownames(scoresX2) <- paste("X2.", 1:nrow(X2), sep = "")
    scores <- rbind(scores, scoresX2)
  }
  
  if(!nrow(plsp$transf$Xscale)){
    plsp$transf$Xscale <- matrix(1, 1,length(plsp$transf$Xcenter))
  }
  
  fresults <- list(scores = scores, 
                   X.loadings = X.loadings, 
                   Y.loadings = Y.loadings, 
                   weights = weights, 
                   projectionM = plsp$projectionM, 
                   variance = list(x.var = exv, y.var = yex), 
                   sc.sdv =cSds(scores), 
                   n.components = nPf, pcSelection = psel,
                   center = plsp$transf$Xcenter,
                   scale = plsp$transf$Xscale)
  fresults$method <- "pls"
  if(pcSel == "opc") {fresults$opcEval <- results} 
  class(fresults) <- c("orthoProjection","list") 
  return(fresults)
}


#' @rdname orthoProjection      
#' @aliases orthoProjection 
#' @aliases plsProjection 
#' @aliases pcProjection 
#' @aliases predict.orthoProjection
#' @export
predict.orthoProjection <- function(object, newdata,...){
  if(missing(newdata))
    return(object$scores)
  else{
    if(length(grep("pca", object$method)) == 0){
      newdata <- sweep(newdata, MARGIN = 2, FUN = "-", STATS = object$center)
      newdata <- sweep(newdata, MARGIN = 2, FUN = "/", STATS = object$scale)
      return(newdata %*% t(object$X.loadings))
    }else{
      predpoj <- projectpls(projectionm = object$projectionM,
                            ncomp = ncol(object$projectionM), 
                            newdata = newdata,
                            scale = TRUE,
                            Xcenter =  object$center,
                            Xscale =  object$scale)
      
      colnames(predpoj) <- paste("pls", 1:ncol(predpoj), sep = "")
      rownames(predpoj) <- rownames(newdata)
      return(predpoj)
    }
  }
}
