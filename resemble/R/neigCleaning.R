#' @title A function for identifying samples that do not belong to any of the neighbourhoods of a given set of samples (neigCleaning)
#' @description
#' This function can be used to identify the samples in a spectral dataset \eqn{Xr} that do not belong to the neighbourhood of any sample in another spectral dataset \eqn{Xu}. 
#' @usage 
#' neigCleaning(Xr, Xu, 
#'              sm = "pc",
#'              pcSelection = list("cumvar", 0.99),
#'              pcMethod = "svd",
#'              Yr = NULL,
#'              ws,
#'              k0,
#'              center = TRUE,
#'              scaled = TRUE,
#'              k.thr, 
#'              k.dist.thr,
#'              k.range,
#'              returnDiss = FALSE, 
#'              cores = 1)
#' @param Xr input (spectral) \code{matrix} (or \code{data.frame}) in which the neighbours of the samples in \code{Xu} shall be searched.
#' @param Xu input (spectral) \code{matrix} (or \code{data.frame}) containing the samples for  which their neighbours will be searched in  \code{Xr}.  
#' @param sm a character string indicating the spectral dissimilarity metric to be used in the selection of the nearest neighbours of each observation for which a prediction is required (see \code{\link{mbl}}). 
#'        Options are: 
#'        \itemize{
#'        \item{\code{"euclid"}: Euclidean dissimilarity.}
#'        \item{\code{"cosine"}: Cosine dissimilarity.}
#'        \item{\code{"sidF"}: Spectral information divergence computed on the spectral variables.}
#'        \item{\code{"sidD"}: Spectral information divergence computed on the density distributions of the spectra.}
#'        \item{\code{"cor"}: Correlation dissimilarity.}
#'        \item{\code{"movcor"}: Moving window correlation dissimilarity.}
#'        \item{\code{"pc"}: Principal components dissimilarity: Mahalanobis dissimilarity computed on the principal components space.}
#'        \item{\code{"loc.pc"}: Dissimilarity estimation based on local principal components.}
#'        \item{\code{"pls"}: Partial least squares dissimilarity: Mahalanobis dissimilarity computed on the partial least squares space.}
#'        \item{\code{"loc.pls"} Dissimilarity estimation based on local partial least squares.}
#'        }
#'        The \code{"pc"} spectral dissimilarity metric is the default. If the \code{"sidD"} is chosen, the default parameters of the \code{sid} function are used however they cab be modified by specifying them as additional arguments in the \code{\link{mbl}} function. This argument can also be set to \code{NULL}, in such a case, a dissimilarity matrix must be specified in the \code{dissimilarityM} argument of the \code{\link{mbl}} function.
#' @param pcSelection if \code{sm = "pc"}, \code{sm = "loc.pc"}, \code{sm = "pls"} or \code{sm = "loc.pls"} a list which specifies the method to be used for identifying the number of principal components to be retained for computing the Mahalanobis distance of each sample in \code{sm = "Xu"} to the centre of \code{sm = "Xr"}. It also specifies the number of components in any of the following cases: \code{sm = "pc"}, \code{sm = "loc.pc"}, \code{sm = "pls"} and \code{sm = "loc.pls"}. This list must contain two objects in the following order: \itemize{
#'        \item{\code{method}:}{the method for selecting the number of components. Possible options are:  \code{"opc"} (optimized pc selection based on Ramirez-Lopez et al. (2013a, 2013b). See the \code{\link{orthoProjection}} function for more details;  \code{"cumvar"} (for selecting the number of principal components based on a given cumulative amount of explained variance); \code{"var"} (for selecting the number of principal components based on a given amount of explained variance); and  \code{"manual"} (for specifying manually the desired number of principal components)}
#'        \item{\code{value}:}{a numerical value that complements the selected method. If \code{"opc"} is chosen, it must be a value indicating the maximal number of principal components to be tested (see Ramirez-Lopez et al., 2013a, 2013b). If \code{"cumvar"} is chosen, it must be a value (higher than 0 and lower than 1) indicating the maximum amount of cumulative variance that the retained components should explain. If \code{"var"} is chosen, it must be a value (higher than 0 and lower than 1) indicating that components that explain (individually) a variance lower than this treshold must be excluded. If \code{"manual"} is chosen, it must be a value specifying the desired number of principal components to retain.
#'        }}
#'        The default method for the \code{pcSelection} argument is \code{"opc"} and the maximal number of principal components to be tested is set to 40.
#'        Optionally, the \code{pcSelection} argument admits \code{"opc"} or \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character string. In such a case the default for \code{"value"} when either \code{"opc"} or \code{"manual"} are used is 40. When \code{"cumvar"} is used the default \code{"value"} is set to 0.99 and when \code{"var"} is used the default \code{"value"} is set to 0.01.
#' @param pcMethod a character string indicating the principal component analysis algorithm to be used. Options are: \code{"svd"} (default) and \code{"nipals"}. See \code{\link{orthoDiss}}.
#' @param Yr either if the method used in the \code{pcSelection} argument is \code{"opc"} or if the \code{sm} argument is either \code{"pls"} or \code{"loc.pls"}, then it must be a \code{vector} containing the side information corresponding to the spectra in \code{Xr}. It is equivalent to the \code{sideInf} parameter of the \code{\link{simEval}} function. It can be a numeric \code{vector} or \code{matrix} (regarding one or more continuous variables). The root mean square of differences (rmsd) is used for assessing the similarity between the samples and their corresponding most similar samples in terms of the side information provided. When \code{sm = "pc"}, this parameter can also be a single discrete variable of class \code{factor}. In such a case the kappa index is used. See \code{\link{simEval}} function for more details.
#' @param ws an odd integer value which specifies the window size when the moving window correlation similarity/dissimilarity is used (i.e \code{sm = "movcor"}). The default value is 41.
#' @param k0 if any of the local similarity/dissimilarity methods is used (i.e. either \code{sm = "loc.pc"} or \code{sm = "loc.pls"}) a numeric integer value. This argument controls the number of initial neighbours(\eqn{k0}) to retain in order to compute the local principal components (at each neighbourhood).
#' @param center a logical indicating if \code{Xr} and \code{Xu} must be centered (on the basis of \eqn{Xr \cup Xu}).
#' @param scaled a logical indicating if \code{Xr} and \code{Xu} must be scaled (on the basis of \eqn{Xr \cup Xu}).
#' @param k.thr an integer value indicating the k-nearest neighbours of each sample in \code{Xu} that must be selected from \code{Xr}. 
#' @param k.dist.thr an integer value indicating a distance treshold. When the distance between a sample in \code{Xr} and a sample in \code{Xu} is below the given treshold, the sample in sample in \code{Xr} is retained, otherwise it is ignored.  The treshold depends
#'        on the corresponding similarity/dissimilarity metric specified in \code{sm}. Either \code{k.thr} or \code{k.dist.thr} must be specified.
#' @param k.range a vector of length 2 which specifies the minimum (first value) and the maximum (second value) number of neighbours allowed when the \code{k.dist.thr} argument is used. 
#' @param returnDiss a logical indicating if the similarity/dissimilarity matrix must be returned. Default is \code{FALSE}.
#' @param cores number of cores used when \code{method} in \code{pcSelection} is \code{"opc"} (which can be computationally intensive) (default = 1).
#' @details
#' This function may be specially useful when the reference set (\code{Xr}) is very large. In some cases the number of observations in the reference set can be reduced by removing irrelevant samples (i.e. samples that are not neighbours of a particular target set). If \code{Xr} is very large, it is recommended to consider the use this function prior using the \code{mbl} function.
#' @return \code{neigCleaning} returns a \code{list} containing the following objects:
#' \itemize{
#'  \item{\code{select}}{ the indices of the observations in \code{Xr} that belong to the negihborhood of the samples in \code{Xu}}.
#'  \item{\code{reject}}{ the indices of the observations in \code{Xr} that do not belong to the negihborhood of the samples in \code{Xu}}.
#'  \item{\code{rn.lower.k.dist}}{ a \code{data.frame} that is returned only if the \code{k.dist.thr} argument was used. It comprises three columns, the first one (\code{sampleIndex}) indicates the index of the samples in \code{Xu}, the second column (\code{nk}) indicates the number of neighbours found in \code{Xr} for each sample in \code{Xr} and the third column (\code{neighbours.used}) indicates whether the original number of neighbours (below the distance treshold) was used or if the number of neighbours was reset to one of the range values specified in the \code{k.range} argument.}
#'  \item{\code{dissimilarity}}{ the distance matrix used.}
#'  }
#' @author Leonardo Ramirez-Lopez
#' @references 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local approach for modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196, 268-279.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{fDiss}}, \code{\link{corDiss}}, \code{\link{sid}}, \code{\link{orthoDiss}}, \code{\link{mbl}}
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
#' # Identify the non-neighbour samples using the default parameters
#' # (In this example all the samples in Xr belong at least to the 
#' # first 100 neighbours of one sample in Xu)
#' ex1 <- neigCleaning(Xr = Xr, Xu = Xu, 
#'                             k.thr = 100)
#' 
#' # Identify the non-neighbour samples using principal component(PC) 
#' # and partial least squares (PLS) distances, and using the "opc" 
#' # approach for selecting the number of components
#' ex2 <- neigCleaning(Xr = Xr, Xu = Xu, 
#'                             Yr = Yr,
#'                             sm = "pc",
#'                             pcSelection = list("opc", 40),
#'                             k.thr = 150)
#' 
#' ex3 <- neigCleaning(Xr = Xr, Xu = Xu, 
#'                             Yr = Yr,
#'                             sm = "pls",
#'                             pcSelection = list("opc", 40),
#'                             k.thr = 150)
#' 
#' # Identify the non-neighbour samples using distances computed 
#' # based on local PC analysis and using the "cumvar" and "var" 
#' # approaches for selecting the number of PCs
#' ex4 <- neigCleaning(Xr = Xr, Xu = Xu, 
#'                             sm = "loc.pc",
#'                             pcSelection = list("cumvar", 0.999),
#'                             k0 = 200,
#'                             k.thr = 150)
#' 
#' ex5 <- neigCleaning(Xr = Xr, Xu = Xu, 
#'                             sm = "loc.pc",
#'                             pcSelection = list("var", 0.001),
#'                             k0 = 200,
#'                             k.thr = 150)
#' }                          
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

neigCleaning <- function(Xr, Xu, 
                         sm = "pc",
                         pcSelection = list("cumvar", 0.99),
                         pcMethod = "svd",
                         Yr = NULL,
                         ws,
                         k0,
                         center = TRUE,
                         scaled = TRUE,
                         k.thr, 
                         k.dist.thr,
                         k.range,
                         returnDiss = FALSE, 
                         cores = 1){

  # Sanity checks
  match.arg(sm, c("euclid", "cosine", "sidF", "sidD", "cor", "movcor", "pc", "loc.pc", "pls", "loc.pls"))
  
  if(!is.numeric(cores))
    stop("The 'cores' argument must be numeric")
  
  if(!is.logical(center))
    stop("'center' argument must be logical")
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical")
  
  if(sm == "movcor"){
    if(missing(ws))
      stop("when 'movcor' is used a moving window size ('ws') must be specified")
    if(ws < 3 | ws > (ncol(Xr) - 1) | length(ws) != 1) 
      stop(paste("ws must be an unique value between 3 and", (ncol(Xr) - 1))) 
    if((ws %% 2) == 0)
      stop("ws must be an odd value")
  }
  
  if(!missing(k.thr) & !missing(k.dist.thr)) 
    stop("Only one of k.thr or k.dist.thr can be specified")  
  
  if(!missing(k.thr)){
    k <- as.integer(k.thr)
    if(k.thr < 10)
      stop("Minimum number of nearest neighbours allowed is 10")
    if(k.thr > nrow(Xr))
      stop("The number of nearest neighbours cannot exceed the number of reference observations")
    if(sm ==  "loc.pc"){
      if(missing(k0))
        stop("For the 'loc.pc' similarity/dissimilarity, the 'k0' argument must be specified")
      if(k0 < k.thr)
        stop("k0 must be larger than 'k.thr'")
    }
  }
  
  if(!missing(k.dist.thr)){
    if(missing(k.range))
      stop("if the k.dist.thr argument is used, k.range must be specified")
    if(length(k.range) != 2 | !is.numeric(k.range) | diff(k.range) < 0)
      stop("k.range must be a vector (of length 2) which specifies the minimum (first value) and the maximum (second value) number of neighbours") 
    k.min <- as.integer(k.range[1])
    k.max <- as.integer(k.range[2])
    if(k.min < 10)
      stop("Minimum number of nearest neighbours allowed is 10")
    if(k.max > nrow(Xr))
      stop("Maximum number of nearest neighbours cannot exceed the number of reference observations")
    if(sm ==  "loc.pc"){
      if(missing(k0))
        stop("For the 'loc.pc' similarity/dissimilarity, the 'k0' argument must be specified")
      if(k0 < k.max)
        stop("k0 must be larger than the maximum number of nearest neighbours in the 'k.range' argument")
    }
  } 
  
  if(sm == "euclid"){ 
    d.mat <- fDiss(Xr = Xr, X2 = Xu, method = "euclid", center = FALSE, scaled = FALSE)
    d.mat <- d.mat # standardize by the number of columns
  } 
  
  if(sm == "pc"){
    mtd <- ifelse(pcMethod == "svd", "pca", "pca.nipals")  
    orthD <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = pcSelection, Yr, method = mtd, center = center, scaled = scaled, return.all = FALSE, call. = FALSE, cores = cores)
    d.mat <- orthD$dissimilarity
  }
  
  if(sm == "loc.pc"){
    mtd <- ifelse(pcMethod == "svd", "pca", "pca.nipals")  
    orthD <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = pcSelection, Yr, method =  mtd, local = TRUE, k0 = k0, center = center, scaled = scaled, return.all = FALSE, call. = FALSE)
    d.mat <- orthD$dissimilarity[ , , drop = FALSE]
  }
  
  if(sm == "pls"){
    orthD <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = pcSelection, Yr, method = "pls", center = center, scaled = scaled, return.all = FALSE, call. = FALSE, cores = cores)
    d.mat <- orthD$dissimilarity
  }
  
  if(sm == "loc.pls"){
    orthD <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = pcSelection, Yr, method = "pls", local = TRUE, k0 = k0, center = center, scaled = scaled, return.all = FALSE, call. = FALSE)
    d.mat <- orthD$dissimilarity[ , , drop = FALSE]
  }
  
  if(sm == "movcor") {
    d.mat <- corDiss(Xr = Xr, X2 = Xu, ws = ws, center = FALSE, scaled = FALSE)
  }
  if(sm == "cor") {
    d.mat <- corDiss(Xr = Xr, X2 = Xu, ws = NULL, center = FALSE, scaled = FALSE)
  }      
  
  if(sm == "sidF") {
    d.mat <- sid(Xr = Xr, X2 = Xu, mode = "feature", center = FALSE, scaled = FALSE, reg = 10^-4)$sid
  }
  if(sm == "sidD") {
    d.mat <- sid(Xr = Xr, X2 = Xu, mode = "density", center = FALSE, scaled = FALSE, reg = 10^-4)$sid
  }
  if(sm == "cosine") {
    d.mat <- fDiss(Xr = Xr, X2 = Xu, method = "cosine", center = FALSE, scaled = FALSE)
  }
  
  ####################################################
  
  dimnames(d.mat) <- list(1:nrow(d.mat), 1:ncol(d.mat))
  
  if(!missing(k.dist.thr)){
    da <- d.mat  
    da[d.mat <= k.dist.thr] <- 1
    da[d.mat > k.dist.thr] <- 0
    
    n.k <- colSums(da)
    
    dda <- da[,!(n.k < k.min | n.k > k.max)]
    
    if((length(dda) !=0)){
      if(is.vector(dda))
        dda <- as.matrix(dda)
      select <- (1:nrow(d.mat))[rowSums(dda) != 0]
    }else {select <- NULL}
    
    if(sum(n.k < k.min) != 0)
    {
      select <- c(select, (apply(as.matrix(d.mat[,n.k < k.min]), 2, order)[1:k.min,]))
    }
    
    if(sum(n.k > k.max) != 0)
    {
      select <- c(select, (apply(as.matrix(d.mat[,n.k > k.max]), 2, order)[1:k.max,]))
    }
    
    select <- unique(select)
    reject <- (1:nrow(d.mat))[!(1:nrow(d.mat) %in% select)]
    
    repl <- rep("original value", length(n.k))
    repl[n.k > k.max] <- paste("reset to", k.max)
    repl[n.k < k.min] <- paste("reset to", k.min)
    
    if(length(reject) == 0)
      reject <- "No samples to reject"
    
    if(length(select) == 0)
      select <- "No samples to select"
    
    rslt <- list(select = sort(select), reject = reject, n.lower.k.dist = data.frame(sampleIndex = 1:length(n.k),nk = n.k, neighbours.used = repl))
    if(returnDiss)
      rslt$dissimilarity <- d.mat
    return(rslt)
  }else{
    select <- sort(unique(as.vector(apply(d.mat, 2, order)[1:k.thr,])))
    reject <- (1:nrow(d.mat))[!(1:nrow(d.mat) %in% select)]
  }
  
  if(length(reject) == 0)
    reject <- "No samples to reject"
  
  if(length(select) == 0)
    select <- "No samples to select"

  rslt <- list(select = sort(select), reject = reject)
  if(returnDiss)
    rslt$dissimilarity <- d.mat
  return(rslt)
}
