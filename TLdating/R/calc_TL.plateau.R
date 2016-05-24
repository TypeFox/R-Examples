#' Estimate value for plateau test
#' 
#' The function estimates the quotient between the natural and the
#' additive/regenerate signal.
#' 
#' 
#' @param Ln
#'  \link{numeric} (\bold{required}): Ln vector
#' @param Ln.error
#'  \link{numeric} (\bold{required}): Error for the Ln vector
#' @param Lx 
#'  \link{numeric} (\bold{required}): Ln matrix
#' @param Lx.error 
#'  \link{numeric} (\bold{required}): Error for the Lx matrix
#' 
#' @return
#'  The function provides an \linkS4class{TLum.Results} object containing: \cr
#'  \describe{
#'    \item{\code{LnLx}}{
#'      \link{numeric}: Ln/Lx matrix}
#'    \item{\code{LnLx.error}}{
#'      \link{numeric}: Error for the Ln/Lx matrix.}
#'  }
#'   
#' @author David Strebler, University of Cologne (Germany).
#' 
## @export calc_TL.plateau

calc_TL.plateau <- function(

  Ln,
  
  Ln.error,
  
  Lx,
  
  Lx.error
  
){
  # Integrity Check ---------------------------------------------------------
  if (missing(Ln)){
    stop("[calc_TL.plateau] Warning: Error is missing.")
  }
  if (missing(Ln.error)){
    stop("[calc_TL.plateau] Warning: Error is missing.")
  }
  if (missing(Lx)){
    stop("[calc_TL.plateau] Warning: Error is missing.")
  }
  if (missing(Lx.error)){
    stop("[calc_TL.plateau] Warning: Error is missing.")
  }

  if (!is(Ln,"numeric")){
    stop("[calc_TL.plateau] Error: Input Ln is not of type 'numeric'.")
  } 
  if (!is(Ln.error,"numeric")){
    stop("[calc_TL.plateau] Error: Input Ln.error is not of type 'numeric'.")
  } 
  if (!is(Lx,"matrix")){
    stop("[calc_TL.plateau] Error: Input Lx is not of type 'matrix'.")
  } 
  if (!is(Lx.error,"matrix")){
    stop("[calc_TL.plateau] Error: Input Lx.error is not of type 'matrix'.")
  }  
  if(length(Ln) != length(Ln.error)){
    stop("[calc_TL.plateau] Error: Ln and Ln.error have a different length.")
  }
  if(length(Lx) != length(Lx.error)){
    stop("[calc_TL.plateau] Error: Lx and Lx.error have a different length.")
  }
  if(!is.null(nrow(Lx)) && length(Ln) != nrow(Lx)){
    stop("[calc_TL.plateau] Error: Ln and Lx have a different length.")
  }
  
  # -------------------------------------------------------------------------------
  # Check
  if(length(Ln) == 0){
    stop("[calc_TL.plateau] Error: No natural signal.")
  }
  
  if(length(Lx) == 0){
    stop("[calc_TL.plateau] Error: No additive/regenerated signal.")
  }
  # -------------------------------------------------------------------------------
  
  Ln.error.r <- abs(Ln.error/Ln)
  
  LnLx <- vector()
  LnLx.error <- vector()
  
  for (i in 1:ncol(Lx)){
    temp.Lx <- Lx[,i]
    temp.error <- Lx.error[,i]
    temp.error.r <- temp.error/temp.Lx
    
    # -------------------------------------------------------------------------------
    # Check
    if(length(Ln) != length(temp.Lx)){
      stop("[calc_TL.plateau] Error: The natural and the additive/regenerated curves have different length.")
    }
    # -------------------------------------------------------------------------------    
    
    temp.LnLx <- Ln/temp.Lx
    temp.LnLx.error.r <- sqrt(Ln.error.r^2 + temp.error.r^2)
    temp.LnLx.error <- abs(temp.LnLx.error.r*temp.LnLx)
    
    LnLx <- cbind(LnLx, temp.LnLx)
    LnLx.error <- cbind(LnLx.error, temp.LnLx.error)
  }
  
  LnLx[!is.finite(LnLx)]<- NA
  LnLx.error[!is.finite(LnLx.error)]<- NA
  
  colnames(LnLx) <- colnames(Lx)
  colnames(LnLx.error) <- colnames(Lx.error)
    
  result <- list(LnLx=LnLx,
                 Error=LnLx.error
                 )
  
  new.TLum.Results.calc_TL.plateau <- set_TLum.Results(data = result)
  
  return (new.TLum.Results.calc_TL.plateau)
}
