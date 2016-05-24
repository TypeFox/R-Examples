#' @title A function for evaluating similarity/dissimilarity matrices (simEval) 
#' @description
#' This function searches for the most similar sample of each sample in a given data set based on 
#' a similarity/dissimilarity (e.g. distance matrix). The samples are compared against their corresponding 
#' most similar samples in terms of the side information provided. The root mean square of differences and 
#' the correlation coefficient are computed for continuous variables and for discrete variables the kappa index 
#' is calculated. 
#' @usage 
#' simEval(d, sideInf, lower.tri = FALSE, cores = 1, ...)
#' @param d a \code{vector} or a square symmetric \code{matrix} (or \code{data.frame}) of similarity/dissimilarity scores between samples of a given dataset (see \code{lower.tri}).
#' @param sideInf a \code{vector} containing the side information corresponding to the samples in the dataset from which the similarity/dissimilarity matrix was computed. It can be either a numeric vector (continuous variable) or a factor (discrete variable). If it is a numeric \code{vector},  the root mean square of differences is used for assessing the similarity between the samples and their corresponding most similar samples in terms of the side information provided. If it is a factor, then the kappa index is used. See details.
#' @param lower.tri a \code{logical} indicating whether the input similarities/dissimilarities are given as a \code{vector} of the lower triangle of the distance matrix (as returned e.g. by \code{base::dist}) or as a square symmetric \code{matrix} (or \code{data.frame}) (default = \code{FALSE})
#' @param cores number of cores used to find the neareast neighbours of similarity/dissimilarity scores (default = 1). See details.
#' @param ... additional parameters (for internal use only).
#' @details
#' For the evaluation of similarity/dissimilarity matrices this function uses side information (information about one variable which is available for a
#' group of samples, Ramirez-Lopez et al., 2013). It is assumed that there is a correlation (or at least an indirect or secondary correlation)
#' between this side informative variable and the spectra. In other words, this approach is based on the assumption that the similarity measures between the spectra of a given group of samples should be able to reflect their
#' similarity also in terms of the side informative variable (e.g. compositional similarity). 
#' If \code{sideInf} is a numeric \code{vector} the root mean square of differences (RMSD) is used for assessing the similarity between the samples and their corresponding most similar samples in terms of the side information provided. It is computed as follows:
#' It can be computed as:
#' \deqn{RMSD = \sqrt{\frac{1}{n} \sum_{i=1}^n {(y_i - \ddot{y}_i)^2}}}
#' where \eqn{y_i} is the value of the side variable of the \eqn{i}th sample, \eqn{\ddot{y}_i} is the value of the side variable of the nearest neighbour 
#' of the \eqn{i}th sample and \eqn{n} is the total number of observations. 
#' If \code{sideInf} is a factor the kappa index (\eqn{\kappa}) is used instead the RMSD. It is computed as follows:
#' \deqn{\kappa = \frac{p_{o}-p_{e}}{1-p_{e}}}
#' where both \eqn{p_o} and \eqn{p_e} are two different agreement indexes between the the side information of the samples and the side information of their corrresponding nearest samples (i.e. most similar samples). 
#' While \eqn{p_o} is the relative agreement \eqn{p_e} is the the agreement expected by chance. 
#' Multi-threading for the computation of dissimilarities (see \code{cores} parameter) is based on OpenMP and hence works only on windows and linux. 
#' @return \code{simEval} returns a \code{list} with the following components:
#' \itemize{
#'  \item{"\code{eval}}{ either the RMSD (and the correlation coefficient) or the kappa index}
#'  \item{\code{firstNN}}{ a \code{data.frame} containing the original side informative variable in the first column and the side informative values of the corresponding nearest neighbours in the second column}
#'  }
#' @author Leonardo Ramirez-Lopez
#' @references 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local approach for modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196, 268-279.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @examples
#' \dontrun{
#' require(prospectr)
#' 
#' data(NIRsoil)
#' 
#' Yr <- NIRsoil$Nt[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' # Example 1
#' # Compute a principal components distance
#' pca.d <- orthoDiss(Xr = Xr, pcSelection = list("cumvar", 0.999), 
#'                    method = "pca", 
#'                    local = FALSE, 
#'                    center = TRUE, scaled = TRUE)
#' 
#' # The final number of pcs used for computing the distance 
#' # matrix of objects in Xr
#' pca.d$n.components
#' 
#' # The final distance matrix 
#' ds <- pca.d$dissimilarity
#' 
#' # Example 1.1
#' # Evaluate the distance matrix on the baisis of the 
#' # side information (Yr) associated with Xr
#' se <- simEval(d = ds, sideInf = Yr)
#' 
#' # The final evaluation results
#' se$eval
#' 
#' # The final values of the side information (Yr) and the values of 
#' # the side information corresponding to the first nearest neighbours 
#' # found by using the distance matrix
#' se$firstNN
#' 
#' # Example 1.2
#' # Evaluate the distance matrix on the baisis of two side 
#' # information (Yr and Yr2) 
#' # variables associated with Xr
#' Yr2 <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' se2 <- simEval(d = ds, sideInf = cbind(Yr, Yr2))
#' 
#' # The final evaluation results
#' se2$eval
#' 
#' # The final values of the side information variables and the values 
#' # of the side information variables corresponding to the first 
#' # nearest neighbours found by using the distance matrix
#' se2$firstNN
#' 
#' ###
#' # Example 2
#' # Evaluate the distances produced by retaining different number of 
#' # principal components (this is the same principle used in the 
#' # optimized principal components approach ("opc"))
#' 
#' # first project the data
#' pca <- orthoProjection(Xr = Xr, method = "pca", 
#'                        pcSelection = list("manual", 30), 
#'                        center = TRUE, scaled = TRUE)
#' 
#' # standardize the scores
#' scores.s <- sweep(pca$scores, MARGIN = 2, 
#'                   STATS = pca$sc.sdv, FUN = "/")
#' rslt <-  matrix(NA, ncol(scores.s), 3)
#' colnames(rslt) <- c("pcs", "rmsd", "r")
#' rslt[,1] <- 1:ncol(scores.s)
#' for(i in 1:ncol(scores.s))
#' {
#'   sc.ipcs <- scores.s[ ,1:i, drop = FALSE]
#'   di <- fDiss(Xr = sc.ipcs, method = "euclid", 
#'               center = FALSE, scaled = FALSE)
#'   se <- simEval(d = di, sideInf = Yr)
#'   rslt[i,2:3] <- unlist(se$eval)
#' }
#' plot(rslt) 
#' 
#' ###
#' # Example 3
#' # Example 3.1
#' # Evaluate a dissimilarity matrix computed using a moving window 
#' # correlation method
#' mwcd <- mcorDiss(Xr = Xr, ws = 35, center = FALSE, scaled = FALSE)
#' se.mw <- simEval(d = mwcd, sideInf = Yr)
#' se.mw$eval
#' 
#' # Example 3.2
#' # Evaluate a dissimilarity matrix computed using the correlation 
#' # method
#' cd <- corDiss(Xr = Xr, center = FALSE, scaled = FALSE)
#' se.nc <- simEval(d = cd, sideInf = Yr)
#' se.nc$eval
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
## 13.03.2014 Antoine The explanation of the cores argument was modified
## 18.03.2014 Antoine Add error message when input dissimilarity matrix 
##                    in simEval is not squared


simEval <- function(d, sideInf, lower.tri = FALSE, cores = 1, ...){
 
  in.call <- match.call()
  if(!is.null(in.call$call.))
    call. <- in.call$call.
  else
    call. <- TRUE
  
  if(!is.numeric(cores))
    stop("The 'cores' argument must be numeric")
  if(!is.logical(lower.tri))
    stop("'lower.tri' argument must be logical")
  if(is.numeric(as.matrix(sideInf))){
    sideInf <- as.matrix(sideInf)
    ny <- ncol(sideInf)
    if(!lower.tri){
      if(!(is.data.frame(d) | is.matrix(d)))
        stop("'d' must be a matrix or a data.frame when lower.tri = FALSE", call. = call.)
      dimd <- dim(d)
      if(dimd[1] != dimd[2])
        stop("'d' must be a square matrix when lower.tri = FALSE", call. = call.)
      if(nrow(d) != nrow(sideInf))
        stop("The number of rows of the 'd' matrix does not match the number of observations in 'sideInf'", call. = call.)
      most <- which_min(d,cores)  # find nearest neighbours
    } else {
      if(!ifelse(is.matrix(d),ncol(d)==1,is.vector(d)|is.atomic(d)))
        stop("'d' must be a vector or a matrix of 1 column when lower.tri = TRUE", call. = call.)
      if(length(d) != (nrow(sideInf)^2-nrow(sideInf))/2)
        stop("The length the 'd' vector does not match the number of observations in 'sideInf'", call. = call.)
      most <- which_minV(d,cores)
    }
    if(sum(colSums(!is.na(sideInf)) < 4) > 0) 
      stop("At least one of the side information variables contains less than 4 values", call. = call.)
    rslt <- data.frame(rmsd = rep(NA, ncol(sideInf)), r = rep(NA, ncol(sideInf)))
    if(!is.null(colnames(sideInf)))
    {
      rownames(rslt) <- colnames(sideInf)
    }else{
      rownames(rslt) <- paste("sideInf.", 1:ncol(sideInf), sep ="")
      colnames(sideInf) <- rownames(rslt)
    }
    for(i in 1:ncol(sideInf))
    {
      rslt$rmsd[i] <- (mean((sideInf[most,i] - sideInf[,i])^2, na.rm = TRUE))^0.5
      rslt$r[i] <- cor(sideInf[most,i][!is.na(sideInf[,i])],sideInf[!is.na(sideInf[,i]),i],use = "pairwise.complete.obs")
    }
    if(ny >1){
      mn.s <- data.frame(mn.sd.rmsd = mean(rslt$rmsd/cSds(sideInf[complete.cases(sideInf),])), mn.r = mean(rslt$r))
      return(list(eval = rslt, global.eval = mn.s, firstNN = data.frame(sideInf, firstNN = sideInf[most,])))
    }else{
      return(list(eval = rslt, firstNN = data.frame(sideInf, firstNN = sideInf[most,])))
    }
  }else{
    if(!is.factor(sideInf))
      stop("sideInf must be either a numeric or a categorical (factor) variable", call. = call.)
    if(!lower.tri){
      if(!(is.data.frame(d)|is.matrix(d)))
        stop("'d' must be a matrix or a data.frame when lower.tri = FALSE", call. = call.)      
      if(nrow(d) != length(sideInf))
        stop("The number of rows of the 'd' matrix does not match the length of 'sideInf'. Note: In the case of categorical variables, only one variable is allowed", call. = call.)
      most <- which_min(d,cores)  # find nearest neighbours      
      } else {
        if(!ifelse(is.matrix(d),ncol(d)==1,is.vector(d)|is.atomic(d)))
          stop("'d' must be a vector or a matrix of 1 column when lower.tri = TRUE", call. = call.)        
        if(length(d) != (nrow(sideInf)^2-nrow(sideInf))/2)
         stop("The length the 'd' vector does not match the number of observations in 'sideInf'. Note: In the case of categorical variables, only one variable is allowed", call. = call.)
        most <- which_minV(d,cores)  # find nearest neighbours
      }
    if(sum(!is.na(sideInf)) < 4)
      stop("Side information contains less than 4 values", call. = call.)
    tab <- as.matrix(table(sideInf[!is.na(sideInf)], sideInf[most][!is.na(sideInf)]))
    total <- sum(tab)
    tab <- tab/total
    p <- rowSums(tab) %*% t(colSums(tab))
    pra <-  sum(diag(tab))
    pre <-  sum(diag(p))
    kappa <- (pra - pre)/(1 - pre)
    return(list(eval = data.frame(kappa), firstNN = data.frame(sideInf, firstNN = sideInf[most])))
  }
} 
