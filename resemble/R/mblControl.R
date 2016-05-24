#' @title A function that controls some aspects of the memory-based learning process in the \code{mbl} function
#' @description
#' This function is used to specify various aspects in the memory-based learning process in the \code{mbl} function 
#' @usage 
#' mblControl(sm = "pc",
#'            pcSelection = list("opc", 40),
#'            pcMethod = "svd",
#'            ws = if(sm == "movcor") 41,
#'            k0,
#'            returnDiss = FALSE,
#'            center = TRUE,
#'            scaled = TRUE,
#'            valMethod = c("NNv", "loc_crossval"),
#'            localOptimization = TRUE,
#'            resampling = 10, 
#'            p = 0.75,
#'            range.pred.lim = TRUE,
#'            progress = TRUE,
#'            cores = 1,            
#'            allowParallel = TRUE)
#' 
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
#'        The \code{"pc"} spectral dissimilarity metric is the default. If the \code{"sidD"} is chosen, the default parameters of the \code{sid} function are used however they cab be modified by specifying them as additional arguments in the \code{\link{mbl}} function.  
#'        
#'        This argument can also be set to \code{"none"}, in such a case, a dissimilarity matrix must be specified in the \code{dissimilarityM} argument of the \code{\link{mbl}} function.
#' @param pcSelection a list which specifies the method to be used for identifying the number of principal components to be retained for computing the Mahalanobis dissimilarity of each sample in \code{sm = "Xu"} to the centre of \code{sm = "Xr"}. It also specifies the number of components in any of the following cases: \code{sm = "pc"}, \code{sm = "loc.pc"}, \code{sm = "pls"} and \code{sm = "loc.pls"}. This list must contain two objects in the following order: \itemize{
#'        \item{\code{method}:}{the method for selecting the number of components. Possible options are:  \code{"opc"} (optimized pc selection based on Ramirez-Lopez et al. (2013a, 2013b). See the \code{\link{orthoProjection}} function for more details;  \code{"cumvar"} (for selecting the number of principal components based on a given cumulative amount of explained variance); \code{"var"} (for selecting the number of principal components based on a given amount of explained variance); and  \code{"manual"} (for specifying manually the desired number of principal components)}
#'        \item{\code{value}:}{a numerical value that complements the selected method. If \code{"opc"} is chosen, it must be a value indicating the maximal number of principal components to be tested (see Ramirez-Lopez et al., 2013a, 2013b). If \code{"cumvar"} is chosen, it must be a value (higher than 0 and lower than 1) indicating the maximum amount of cumulative variance that the retained components should explain. If \code{"var"} is chosen, it must be a value (higher than 0 and lower than 1) indicating that components that explain (individually) a variance lower than this threshold must be excluded. If \code{"manual"} is chosen, it must be a value specifying the desired number of principal components to retain.
#'        }}
#'        The default method for the \code{pcSelection} argument is \code{"opc"} and the maximal number of principal components to be tested is set to 40.
#'        Optionally, the \code{pcSelection} argument admits \code{"opc"} or \code{"cumvar"} or \code{"var"} or \code{"manual"} as a single character string. In such a case the default for \code{"value"} when either \code{"opc"} or \code{"manual"} are used is 40. When \code{"cumvar"} is used the default \code{"value"} is set to 0.99 and when \code{"var"} is used the default \code{"value"} is set to 0.01.
#' @param pcMethod a character string indicating the principal component analysis algorithm to be used. Options are: \code{"svd"} (default) and \code{"nipals"}. See \code{\link{orthoDiss}}.
#' @param ws an odd integer value which specifies the window size when the moving window correlation dissimilarity is used (i.e \code{sm = "movcor"}). The default is 41.
#' @param k0 if any of the local dissimilarity methods is used (i.e. either \code{sm = "loc.pc"} or \code{sm = "loc.pls"}) a numeric integer value. This argument controls the number of initial neighbours(\eqn{k0}) to retain in order to compute the local principal components (at each neighbourhood).
#' @param returnDiss a logical indicating if the dissimilarity matrices must be returned.
#' @param center a logical indicating whether or not the predictor variables must be centered at each local segment (before regression).
#' @param scaled a logical indicating whether or not the predictor variables must be scaled at each local segment (before regression).
#' @param valMethod a character vector which indicates the (internal) validation method(s) to be used for assessing the global performance of the local models. Possible
#'        options are: \code{"NNv"} and \code{"loc_crossval"}. Alternatively \code{"none"} can be used when corss-validation is not required (see details below).
#' @param localOptimization a logical. If \code{valMethod = "loc_crossval"}, it optmizes the parameters of the local pls models (i.e. pls factors for \code{pls} and minimum and maximum pls factors for \code{wapls1}).
#' @param resampling a value indicating the number of resampling iterations at each local segment when \code{"loc_crossval"} is selected in the \code{valMethod} argument. Default is 10.
#' @param p a value indicating the percentage of samples to be retained in each resampling iteration at each local segment when \code{"loc_crossval"} is selected in the \code{valMethod} argument. Default is 0.75 (i.e. 75 "\%")
#' @param range.pred.lim a logical value. It indicates whether the prediction limits at each local regression are determined by the range of the response variable values employed at each local regression. If \code{FALSE}, no prediction limits are imposed. Default is \code{TRUE}.
#' @param progress a logical indicating whether or not to print a progress bar for each sample to be predicted. Default is \code{TRUE}. Note: In case multicore processing is used, this progress bar will not be printed.
#' @param cores number of cores used for the computation of dissimilarities when \code{method} in \code{pcSelection} is \code{"opc"} (which can be computationally intensive) (default = 1). See details.
#' @param allowParallel To allow parallel execution of the sample loop (default is \code{TRUE})
#' @details
#' The validation methods avaliable for assessing the predictive performance of the memory-based learning method used are described as follows:
#'  \itemize{
#'  \item{Leave-nearest-neighbour-out cross validation (\code{"NNv"}):}{ From the group of neighbours of each sample to be predicted, the nearest sample (i.e. the most similar sample) is excluded and then a local model is fitted using the remaining neighbours. This model is then used to predict the value of the target response variable of the nearest sample. These predicted values are finally cross validated with the actual values (See Ramirez-Lopez et al. (2013a) for additional details). This method is faster than \code{"loc_crossval"}}
#'  \item{Local leave-group-out cross validation (\code{"loc_crossval"}):}{ The group of neighbours of each sample to be predicted is partitioned into different equal size subsets. Each partition is selected based on a stratified random sampling which takes into account the values of the response variable of the corresponding set of neighbours. The selected local subset is used as local validation subset and the remaining samples are used for fitting a model. This model is used to predict the target response variable values of the local validation subset and the local root mean square error is computed. This process is repeated \eqn{m} times and the final local error is computed as the average of the local root mean square error of all the \eqn{m} iterations. In the \code{mbl} function \eqn{m} is controlled by the \code{resampling} argument and the size of the subsets is controlled by the \code{p} argument which indicates the percentage of samples to be selected from the subset of nearest neighbours. The global error of the predictions is computed as the average of the local root mean square errors.}
#'  \item{No validation (\code{"none"}):}{ No validation is carried out. If \code{"none"} is seleceted along with \code{"NNv"} and/or \code{"loc_crossval"}, then it will be ignored and the respective validation(s) will be carried out.} 
#'  }
#' Multi-threading for the computation of dissimilarities is based on OpenMP and hence works only on windows and linux. 
#' However, the loop used to iterate over the \code{Xu} samples in \code{mbl} uses the \code{\%dopar\%} operator of the \code{\link{foreach}} package, which can be used to parallelize this internal loop. The last example given in the \code{\link{mbl}} function ilustrates how to parallelize the \code{\link{mbl}} function.
#' @return \code{mblControl} returns a \code{list} of class \code{mbl} with the specified parameters
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @references 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local approach for modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196, 268-279.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#' @seealso \code{\link{fDiss}}, \code{\link{corDiss}}, \code{\link{sid}}, \code{\link{orthoDiss}}, \code{\link{mbl}}
#' @examples
#' #A control list with the default parameters
#' mblControl()
#' 
#' #A control list which specifies the moving correlation 
#' #dissimilarity metric with a moving window of 30
#' mblControl(sm = "movcor", ws = 31)
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
## 18.11.2015 Leo     A sanity check for avoiding potential sm = NULL was 
##                    introduced
## 15.12.2015 Leo     A bug when checking the valMethod provided was fixed
## 03.01.2016 Leo     the localOptimization argument was added


mblControl <- function(sm = "pc",
                       pcSelection = list("opc", 40),
                       pcMethod = "svd",
                       ws = if(sm == "movcor") 41,
                       k0,
                       returnDiss = FALSE,
                       center = TRUE,
                       scaled = TRUE,
                       valMethod = c("NNv", "loc_crossval"),
                       localOptimization = TRUE,
                       resampling = 10, 
                       p = 0.75,
                       range.pred.lim = TRUE,
                       progress = TRUE, cores = 1,
                       allowParallel = TRUE){
  # Sanity checks
  if(!is.logical(allowParallel))
    stop("allowParallel must be a logical value")
  
  
  ## "none" must go first in case sm = NULL (argument wrongly set)
  sm <- match.arg(sm, c("none", "euclid", "cosine", "sidF", "sidD", "cor", "movcor", "pc", "loc.pc", "pls", "loc.pls"))
  
  if(sm == "none"){
    message("A spectral dissimilarity metric ('sm' argument) was not specified and it has been set to 'none'. You will have to provide a dissimilarity matrix to the mbl function by using the 'dissimilarityM' argument")
  }
  
  
  if(sm == "pc")
    match.arg(pcMethod, c("svd", "nipals"))
  
  if(sm == "movcor"){
    if(ws < 3 | length(ws) != 1) 
      stop(paste("'ws' must be an unique odd value greater than 2")) 
    if((ws %% 2) == 0)
      stop("'ws' must be an odd value")
    smParam <- ws
  }
  
  if(sm %in% c("loc.pc", "loc.pls")){
    if(missing(k0))
      stop("If 'sm' is either 'loc.pc' or 'loc.pls', k0 must be provided")
    if(k0 < 10 | length(k0) != 1) 
      stop(paste("'k0' must be an positive integer value greater than 10")) 
    smParam <- k0
  }
  
  if(!is.logical(returnDiss))
    stop("'returnDiss' must be logical")
  
  pcSel <- match.arg(pcSelection[[1]], c("opc", "var", "cumvar", "manual"))
  
  if(length(pcSelection) > 1){
    if(!is.list(pcSelection))
      stop("The pcSelection argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifying only the method can be used, in this case the parameter value is set automatically")     
    pcSelection <- list(method = pcSelection[[1]], value = pcSelection[[2]])
    names(pcSelection) <- c("method", "value")
    if(!is.numeric(pcSelection[[2]])){
      stop("The second object in 'pcSelection' must be a numeric value")
    }
  }
  
  if(pcSel == "opc")
  {
    if(is.list(pcSelection)){
      if(pcSelection$value < 2) 
        stop(paste("The number of principal components must be an integer value greater than 1")) 
    } else{
      pcSelection <- list(method = "manual", value = 40)
      message(paste("The number of principal components to be tested in the 'opc' method was set to 40.", "\n", "Note: An user-defined number of principal components to be retained can be specified in the pcSelection argument with a list in which the first object indicates the method 'manual' and the second object indicates the number of principal components"))
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
        stop("The 'pcSelection' argument must be a list in which the first object indicates the selection method and the second object indicates the parameter value of the method. Optionally, instead a list, a character string specifying only the method can be used, in this case the parameter value is set automatically")     
      pcSelection <- list(pcSelection[[1]], pcSelection[[2]])
      names(pcSelection) <- c("method", "value")
      if(!is.numeric(pcSelection$value)) 
        stop("The second object in 'pcSelection' must be a numeric value")
      if(pcSelection$value > 1 | pcSelection$value <= 0) 
        stop(paste("When the method for 'pcSelection' is either 'var' or 'cumvar' the value in 'pcSelection' must be a number higher than 0 and lower than/or equal to 1"))
    }
  }
  
  if(pcSel == "manual")
  {
    if(is.list(pcSelection)){
      if(pcSelection$value < 2) 
        stop(paste("the number of principal components must be an integer value greater than 1")) 
    } else{
      pcSelection <- list(method = "manual", value = 3)
      message(paste("the number of principal components to be retained was set to 3.", "\n", "Note: An user-defined number of principal components to be retained can be specified in the pcSelection argument with a list in which the first object indicates the method 'manual' and the second object indicates the number of principal components"))
    }
  }
  
  if(!is.logical(center))
    stop("'center' argument must be logical")
  
  if(!is.logical(scaled))
    stop("'scaled' argument must be logical")
  
  if(sum(valMethod %in% c("NNv", "loc_crossval", "none")) != length(valMethod))
    stop("'valmethod' must be one at least one of 'NNv', 'loc_crossval', 'none'")
  
  if("none" %in% valMethod)
    if(length(valMethod) > 1)
      valMethod <- valMethod[!(valMethod == "none")]
  
  if("loc_crossval" %in% valMethod)
  {
    if(!is.numeric(resampling) | length(resampling) != 1)
      stop("The 'resampling' argument must be a single numeric value")
    
    if(!is.numeric(p) | length(p) != 1 | p >= 1 | p <= 0)
      stop("p must be a single numeric value higher than 0 and lower than 1")
  }
  
  if(!is.logical(range.pred.lim))
    stop("'range.pred.lim' must be logical")
  
  if(!is.logical(progress))
    stop("The 'progress' argument must be logical")
  
  if(!is.numeric(cores))
    stop("The 'cores' argument must be numeric")  
  
  if(sm %in% c("movcor", "loc.pc", "loc.pls")){
    cntrl <- list(sm = sm,
                  smParam = smParam,
                  returnDiss = returnDiss,
                  pcSelection = pcSelection,
                  pcMethod = pcMethod,
                  center = center,
                  scaled = scaled,
                  valMethod = valMethod,
                  localOptimization = localOptimization,
                  resampling = resampling, 
                  p = p,
                  range.pred.lim = range.pred.lim,
                  progress = progress, 
                  cores = cores, allowParallel = allowParallel)
    nm <- ifelse(sm == "movcor", "ws", "k0")
    names(cntrl)[names(cntrl) == "smParam"] <- nm
    
  } else{
    cntrl <- list(sm = sm,
                  returnDiss = returnDiss,
                  pcSelection = pcSelection,
                  pcMethod = pcMethod,
                  center = center,
                  scaled = scaled,
                  valMethod = valMethod,
                  localOptimization = localOptimization,
                  resampling = resampling, 
                  p = p,
                  range.pred.lim = range.pred.lim,
                  progress = progress, 
                  cores = cores, allowParallel = allowParallel)
  }
  return(cntrl)
}
