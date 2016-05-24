#' Separate the additive and the regenerative curves 
#' 
#' Internal function called by \link{analyse_TL.MAAD}. \cr
#' This function separates the additive curves from the regenerative curves using a vector containing the data type associate with each curve. 
#' Additive curves have "Natural" or "N+dose" as datatype.
#' Regenerative curves have "Bleach" or "Bleach+dose" as datatype.
#' Other datatypes are not supported.   
#' 
#'  
#' @param doses
#'  \link{numeric} (\bold{required}): doses vector
#' @param Lx
#'  \link{numeric} (\bold{required}): Lx matrix
#' @param Lx.error
#'  \link{numeric} (\bold{required}): Error for the Lx matrix.
#' @param dTypes
#'  \link{character} (\bold{required}): data type vector.
#' 
#' @author David Strebler, University of Cologne (Germany).
#' 
## @export calc_TL.MAAD.separate

calc_TL.MAAD.separate <- function(
  Lx,
  
  Lx.error,
  
  doses,
  
  dTypes
  
){
  # ------------------------------------------------------------------------------  
  # Integrity Check
  # ------------------------------------------------------------------------------  
  
  if(!is.character(dTypes)){
    stop("[calc_TL.MAAD.separate] Error: names is not of type 'character'.")
  }
  
  if(!is.numeric(doses)){
    stop("[calc_TL.MAAD.separate] Error: Lx is not of type 'numeric'.")
  }
  
  if(!is.numeric(Lx)){
    stop("[calc_TL.MAAD.separate] Error: Lx is not of type 'numeric'.")
  }
  
  if(!is.numeric(Lx.error)){
    stop("[calc_TL.MAAD.separate] Error: Lx.error is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------  
  
  V_ADDITIVE <- c("Natural", "N+dose")
  C_ADDITIVE <- "A"
  
  V_REGENERATIVE <- c("Bleach", "Bleach+dose")
  C_REGENERATIVE <- "R"
  
  # ------------------------------------------------------------------------------  
  #Values Check
  
  if(length(dTypes) != length(doses)){
    stop("[calc_TL.MAAD.separate] Error: dTypes and doses do not have the same size.") 
  }
  
  if(length(Lx) != length(Lx.error)){
    stop("[calc_TL.MAAD.separate] Error: Lx and Lx.error do not have the same size.") 
  }
  
  if(length(doses) != ncol(Lx)){
    stop("[calc_TL.MAAD.separate] Error: Lx and doses do not have the same size.") 
  }
  
  if(FALSE %in% (dTypes %in% c(V_ADDITIVE,V_REGENERATIVE))){
    stop("[calc_TL.MAAD.separate] Error: Data type not supported.") 
  }
  # ------------------------------------------------------------------------------  
  
  uDoses <- unique(doses[order(doses)])
  
  temp.n <- 0
  
  aDoses <- vector()
  aNames <- vector()
  
  aLx <- vector()
  aLx.error <- vector()
  
  rDoses <- vector()
  rNames <- vector()
  
  rLx <- vector()
  rLx.error <- vector()
  
  for(temp.dose in uDoses){
    
    if(temp.dose > 0){
      temp.n <- temp.n+1
    }
    
    #Signal selection based on the dose step
    for(i in 1: length(doses)){
      if(doses[i] == temp.dose){
        if(dTypes[i] %in% V_ADDITIVE){        
          
          temp.Lx <- Lx[,i]
          temp.Lx.error <- Lx.error[,i]
          
          temp.name <- paste(C_ADDITIVE,temp.n, sep="")   
          
          if(temp.name == "A0"){
            temp.name <- "N"
          }
          
          aNames <- c(aNames, temp.name)
          aDoses <- c(aDoses, temp.dose)
          aLx <- cbind(aLx, temp.Lx)
          aLx.error <- cbind(aLx.error, temp.Lx.error)
        
        }else if(dTypes[i] %in% V_REGENERATIVE){
          temp.Lx <- Lx[,i]
          temp.Lx.error <- Lx.error[,i]
          
          temp.name <- paste(C_REGENERATIVE,temp.n, sep="")   
          
          rNames <- c(rNames, temp.name)
          rDoses <- c(rDoses, temp.dose)
          rLx <- cbind(rLx, temp.Lx)
          rLx.error <- cbind(rLx.error, temp.Lx.error)
        }else{
          stop("[calc_TL.MAAD.separate] Error: Data type not supported.")   
        }
      }    
    }
  }  
  
  
  # Column naming
  if(length(aLx)>0){
    colnames(aLx) <- aNames
    colnames(aLx.error) <- aNames    
  }

  if(length(rLx)>0){
    colnames(rLx) <- rNames
    colnames(rLx.error) <- rNames    
  }
  
  #Values Check
  aLx[!is.finite(aLx)] <- NA
  aLx.error[!is.finite(aLx.error)] <- NA
  
  rLx[!is.finite(rLx)] <- NA
  rLx.error[!is.finite(rLx.error)] <- NA
  
  # ---------------------------------------------------------------------
  
  result <- list(aDoses=aDoses,
                 aNames=aNames,
                 aLx=aLx,
                 aLx.error=aLx.error,
                 rDoses=rDoses,
                 rNames=rNames,
                 rLx=rLx,
                 rLx.error=rLx.error)
  
  new.TLum.Results.calc_TL.MAAD.separate <- set_TLum.Results(data = result)
  
  return(new.TLum.Results.calc_TL.MAAD.separate)
}
