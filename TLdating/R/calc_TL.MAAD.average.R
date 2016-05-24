#' Estimate average additive curves for the MAAD protocol.
#' 
#' Internal function called by \link{analyse_TL.MAAD}. \cr
#' This function estimates the average curves for each additive dose.
#' 
#' 
#' @param names
#'  \link{character} (\bold{required}): Names of the TL curves
#' @param doses
#'  \link{numeric} (\bold{required}): additive doses used for the TL curve
#' @param Lx
#'  \link{numeric} (\bold{required}): Lx matrix
#' @param Lx.error  
#'  \link{numeric} (\bold{required}): Error for the Lx matrix.
#' 
#' @return
#'  The function provide an \linkS4class{TLum.Results} object containing: \cr
#'  \describe{
#'    \item{\code{names}}{
#'      \link{character}: Vector with the names of the average additive curves.}
#'    \item{\code{doses}}{
#'      \link{character}: Vector with the additive doses corresponding to each average additive curve.}
#'    \item{\code{Lx}}{
#'      \link{numeric}: new average additive curve Lx matrix.}
#'    \item{\code{Lx.error}}{
#'      \link{numeric}: Error on the new Lx matrix.}
#'  }
#'  
#' @author David Strebler, University of Cologne (Germany).
#' 
## @export calc_TL.MAAD.average

calc_TL.MAAD.average <- function(
  names,
  doses,
  Lx,
  Lx.error  
  
){
  # -----------------------------------
  # Integrity Check
  
  if(missing(names)){
    stop("[calc_TL.MAAD.average] Error: Input object is missing.")
    
  }else if(!is.character(names)){
    stop("[calc_TL.MAAD.average] Error: names is not of type 'character'.")
  }
  
  if(missing(doses)){
    stop("[calc_TL.MAAD.average] Error: Input object is missing.")
    
  }else if(!is.numeric(doses)){
    stop("[calc_TL.MAAD.average] Error: doses is not of type 'numeric'.")
  }
  
  if(missing(Lx)){
    stop("[calc_TL.MAAD.average] Error: Input object is missing.")
    
  }else if(!is.numeric(Lx)){
    stop("[calc_TL.MAAD.average] Error: Lx is not of type 'numeric'.")
  }
  
  if(missing(Lx.error)){
    stop("[calc_TL.MAAD.average] Error: Input object is missing.")
    
  }else if(!is.numeric(Lx.error)){
    stop("[calc_TL.MAAD.average] Error: Lx.error is not of type 'numeric'.")
  }
  
  #------------------------------------
  #Values Check
  
  if(length(names) != length(doses)){
    stop("[calc_TL.MAAD.average] Error: names and doses do not have the same size.") 
  }
  
  if(length(Lx) != length(Lx.error)){
    stop("[calc_TL.MAAD.average] Error: Lx and Lx.error do not have the same size.") 
  }
  
  if(length(doses) != ncol(Lx)){
    stop("[calc_TL.MAAD.average] Error: Lx and doses do not have the same size.") 
  }
  #------------------------------------
  
  new.names <- unique(names)
  new.doses <- unique(doses)
  
  new.Lx <- vector()
  new.Lx.error <- vector()
  
  for(temp.dose in new.doses){
    
    temp.Lx <- vector()
    temp.Lx.error <- vector()
    temp.names <- vector()
    
    #Signal selection based on the dose step
    for(i in 1: length(doses)){
      if(doses[i] == temp.dose){        
        temp.Lx <-cbind(temp.Lx, Lx[,i])
        temp.Lx.error <-cbind(temp.Lx.error, Lx.error[,i])
      }
    }
    
    #weighted average for each dose step      
    temp.w <- 1/(temp.Lx.error^2)
    
    temp.Lx.a <- vector()
    temp.Lx.a.error <- vector()
    
    for(j in 1:nrow(Lx)){
      temp.Lx.a[j] <- sum(temp.w[j,]*temp.Lx[j,],na.rm=TRUE)/sum(temp.w[j,],na.rm=TRUE)
      temp.Lx.a.error[j] <- 1/sqrt(sum(temp.w[j,], na.rm=TRUE))
    }
    
    new.Lx <- cbind(new.Lx, temp.Lx.a)
    new.Lx.error <- cbind(new.Lx.error, temp.Lx.a.error)  
  }
    
  # Column naming
  if(length(new.Lx)>0){
    colnames(new.Lx) <- new.names 
    colnames(new.Lx.error) <- new.names     
  }
  
  #Check values
  new.Lx[!is.finite(new.Lx)] <- NA
  new.Lx.error[!is.finite(new.Lx.error)] <- NA
  # ---------------------------------------------------------------------
  
  result <- list(names=new.names,
                  doses=new.doses,
                  Lx=new.Lx,
                  Lx.error=new.Lx.error)
  
  new.TLum.Results.calc_TL.MAAD.average <- set_TLum.Results(data = result)
  
  return(new.TLum.Results.calc_TL.MAAD.average)
}
