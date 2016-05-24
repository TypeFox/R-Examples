# Calculate NPDE from PDE data
# roxygen comments
#' Calculates individual normalized prediction distribution errors (NPDE) from
#' PDE data.
#'
#' \pkg{nca.npde} calculates individual normalized prediction distribution 
#' errors (NPDE)  of selected NCA metrics from the PDE data.
#' 
#' \pkg{nca.npde} calculates individual normalized prediction distribution 
#' errors (NPDE) of selected NCA metrics from PDE data. The The deviation of
#' each estimated NCA metrics is scaled by the "spread" of the simulated values.
#' By default, this function calculates the NPDE values of all columns of the
#' input data frame.
#' 
#' @param pdedata A data frame containing the prediction distribution errors
#'   (PDE) of NCA metrics
#' @param pdecol The range of column numbers in the data frame containing the
#'   PDE values, which will be used to calculate the corresponding NPDE
#'
#' @return returns the data frame with the NPDE values based on the input data.
#' @export
#'

nca.npde <- function(pdedata,pdecol){
  "qnorm" <- NULL
  rm(list=c("qnorm"))
  
  if (is.null(pdecol)){
    pdecol <- names(pdedata)
  }else if (!is.null(pdecol) && !is.numeric(pdecol)){
    if (!all(pdecol%in%names(pdedata))) stop("All column names given as PDE data column must be present in pdedata")
  }else if (!is.null(pdecol) && is.numeric(pdecol)){
    if (any(pdecol <= 0) | (max(pdecol) > ncol(pdedata))) stop("Column number for PDE data out of range of pdedata.")
    pdecol <- names(pdedata)[pdecol]
  }
  npde                <- pdedata
  npde[,pdecol]       <- lapply(npde[,pdecol], FUN=function(x) qnorm(as.numeric(x)))
  npde[mapply(is.infinite, npde)] <- 0
  
  names(npde)[which(names(npde)%in%pdecol)] <- paste("npde",names(npde)[which(names(npde)%in%pdecol)],sep="")
  return(npde)
}

