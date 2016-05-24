#######
### Function to obtain slices of the synthetic log likelihood
#######
#' Plot slices of the synthetic log-likelihood.
#' 
#' @param object  \code{synlik} object.
#' @param ranges   ranges of values along which we want the slices. If \code{length(parName) == 1} than range has a vector, while
#'                if \code{length(parName) == 2} it have to be a named list of 2 vectors (ex: \code{list("alpha" = 1:10, "beta" = 10:1)}).
#' @param nsim    Number of simulations used to evaluate the synthetic likelihood at each location.
#' @param param  Named vector containing the value of the ALL parameters (including the sliced one). Parameters that are not
#'                in \code{parName} will be fixed to the values in \code{param}.
#' @param pairs if \code{TRUE} the function will produce a 2D slice for every pair of parameters in \code{ranges}. \code{FALSE}
#'              by default.
#' @param draw    If \code{TRUE} the slice will be plotted.
#' @param trans Named vector or list of transformations to be applied to the parameters in \code{parName} 
#'               before plotting {ex: \code{trans = c(s = "exp", d = "exp")}}/ 
#' @param multicore  If \code{TRUE} the \code{object@@simulator} and \code{object@@summaries} functions will
#'                    be executed in parallel. That is the \code{nsim} simulations will be divided in multiple cores.
#' @param ncores  Number of cores to use if \code{multicore == TRUE}.
#' @param cluster An object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be passed to \code{slik()}, see \code{\link{slik}}.
#' @return Either a vector or matrix of log-synthetic likelihood estimates, depending on whether \code{length(parNames) ==} 1 or 2.
#'         These are returned invisibly.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> 
#' @examples
#' data(ricker_sl)
#' 
#'# Plotting slices of the logLikelihood
#'slice(object = ricker_sl, 
#'      ranges = list("logR" = seq(3.5, 3.9, by = 0.01),
#'                    "logPhi" = seq(2, 2.6, by = 0.01),
#'                     "logSigma" = seq(-2, -0.5, by = 0.01)), 
#'      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
#'      nsim = 500)
#'             
#'\dontrun{
#'# Plotting a contour of the logLikelihood
#'slice(object = ricker_sl, 
#'      ranges = list("logR" = seq(3.5, 3.9, by = 0.01),
#'                    "logPhi" = seq(2, 2.6, by = 0.01),
#'                    "logSigma" = seq(-2, -0.5, by = 0.04)), 
#'      pairs = TRUE,
#'      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
#'      nsim = 500, multicore = TRUE)   
#'}                     
#' @export 

slice <- function(object, ranges, nsim, param = object@param, pairs = FALSE, draw = TRUE, trans = NULL, 
                         multicore = FALSE, ncores = detectCores() - 1, cluster = NULL, ...)
{
  if( !is(object, "synlik") ) stop("object has to be of class \"synlik\" ")
  if( length(object@data) == 0 ) stop("There is no data in the object: length(object@data) == 0")
  if( is.null(names(param)) ) names(param) <- names(object@param)
  if( !is.list(ranges) ) stop("ranges has to be a list of named vectors")

  if( !is.null(trans) ) trans <- as.list(trans)
  
  nPar <- length(ranges)
  out <- list()
  
  # Function that will be used by sapply() or clusterApply to evaluate the likelihood
  funToApply <- function(parVal, ...)
  {
    param[parName] <- parVal
    slik(object, param, nsim, multicore = FALSE, cluster = NULL, ...)
  }
  
  if(multicore) {
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "synlik", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated 
    
    # I put the environment of the function to apply to .Global to avoid that exports all the enviroment at every round of lapply.
    environment(funToApply) <- .GlobalEnv
  }
  
  panelDim <- min(ceiling(sqrt(nPar)), 3)
  par(mfrow = c(panelDim, panelDim))
  
  outlik <- list()
  
  if( !pairs )
  {
    # One dimensional slice
    counter <- 1
    for(parName in names(ranges))
    {
      tmpRange <- ranges[[parName]]
      
      if( !is.vector(tmpRange) ) 
        stop( paste("Range of", parName, "has to be a vector") )
      
      if(multicore)
      {
        clusterExport(cluster, "parName", envir = environment())
        llkTransect <- clusterApply(cluster, tmpRange, funToApply, ...)
        
        llkTransect <- unlist(llkTransect)
      } else {
        llkTransect <- sapply(tmpRange, funToApply, ...)
      }
      
      outlik <- c(outlik, llkTransect)
      
      if(draw)
      {
        if( !is.null(trans[[parName]]) )
        {
          tmpRange <- get(trans[[parName]])(tmpRange)
        }
        
        plot(tmpRange, llkTransect, type = 'l', 
             main = parName, 
             xlab = parName, ylab = "Synthetic Log-likelihood")
        
        if( !(counter %% (panelDim^2)) && (counter != nPar) ){ 
          readline(prompt = "Press <Enter> to see the next plot...") 
          counter <- 0
        }
        counter <- counter + 1
      }
    }
  }
  
  # Two dimensional slice
  if( pairs )
  {
    if(length(ranges) < 2) stop("You need at least 2 parameters if pairs == TRUE")
    
    # Create combinations of pairs
    com <- t( .combn(nPar, 2) )
    
    # Calculate likelihood for each pair
    for(ii in 1:nrow(com))
    {
      parName <- names(ranges)[com[ii, ]]
      
      theGrid <- do.call( cbind, expand.grid(ranges[[parName[1]]], ranges[[parName[2]]]) )
      theGrid <- split(t(theGrid), rep(1:nrow(theGrid), each = ncol(theGrid)) )
      nCol <- length(ranges[[parName[2]]])
      nRow <- length(ranges[[parName[1]]])
      llkTransect <- matrix(NA, nRow, nCol)
      
      if(multicore)
      { 
        clusterExport(cluster, "parName", envir = environment())
        tmplik <- clusterApply(cluster, theGrid, funToApply, ...)
      } else {
        tmplik <- sapply(theGrid, funToApply, ...)
      }
      
      kk <- 1
      for(iCol in 1:nCol)
        for(iRow in 0:(nRow-1)) 
        {
          llkTransect[nRow - iRow, iCol] <- tmplik[[kk]]
          kk <- kk + 1
        }
      
      outlik <- c(outlik, llkTransect)
      
      if(draw) 
      {
        if(!is.null(trans))
        {
          ranges[[parName[1]]] <- get(trans[[parName[1]]])(ranges[[parName[1]]])
          ranges[[parName[2]]] <- get(trans[[parName[2]]])(ranges[[parName[2]]])
        }
        
        .plotMatrix(x = llkTransect, title = paste("Transect wrt params", parName[1], "and", parName[2]),
                   xLabels = round(ranges[[parName[2]]], 3), yLabels = rev(round(ranges[[parName[1]]], 3)), xlab = parName[2], 
                   ylab = parName[1], scaleLab = "Log-lik")
      }
    }
  }
  
  if(multicore && clusterCreated) stopCluster(cluster)
  
  return( invisible(llkTransect) )
}
