#' Convert RLum.Analysis object to an TLum.Analysis.
#'
#' This function convert a \linkS4class{RLum.Analysis} object, created using the package'Luminescence', into a \linkS4class{TLum.Analysis} object, which will be used by this package.
#' The TL package is implemented to only works with its own class of object (TLum.Analysis, TLum.Analysis and TLum.Data.Curve).
#'
#' @param object
#'  \code{\linkS4class{RLum.Analysis}} (\bold{required}): object containing the TL curves used for the ED estimation.
#' @param relative.error
#'  \link{numeric} (\bold{required}): Relative error of the TL signals. Generally, it is between 0.02 and 0.1.
#'
#' @details
#' This function use the data from the RLum.Analysis and the relative.error specified to create a absolute error matrix.
#' Then it create a new TLum.Analysis including all the information from the RLum.Analysis and the new absolute error matrix.
#' For practical reason, the TLdating package considers the error as random. It means that the systematic component of the error will be ignored.
#'
#'
#' @author David Strebler, University of Cologne (Germany), \cr David Strebler
#'
#' @export RLum.Analysis2TLum.Analysis


RLum.Analysis2TLum.Analysis <- function(
  object,

  relative.error

){

  # Integrity Check ---------------------------------------------------------
  if(missing(object)){
    stop("[RLum.Analysis2TLum.Analysis] Error: Input object is missing.")

  }else if(!is(object,"RLum.Analysis")){
    stop("[RLum.Analysis2TLum.Analysis] Error: Input object is not of type 'RLum.Analysis'.")
  }

  if(missing(relative.error)){
    stop("[RLum.Analysis2TLum.Analysis] Error: Input relative error is missing.")

  }else if(!is.numeric(relative.error)){
    stop("[RLum.Analysis2TLum.Analysis] Error: Relative error is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Value check
  if(relative.error > 1){
    warning("[RLum.Analysis2TLum.Analysis] Warning: Input 'relative.error' > 1.")

  }else if(relative.error < -1){
    relative.error <- abs(relative.error)
    warning("[RLum.Analysis2TLum.Analysis] Warning: Input 'relative.error' < -1.")

  }else if(relative.error < 0){
    relative.error <- abs(relative.error)
    warning("[RLum.Analysis2TLum.Analysis] Warning: Input 'relative.error' < 0.")
  }
  # ------------------------------------------------------------------------------


  nRecords <- length(object@records)

  new.protocol <- object@protocol
  records <- object@records

  new.records <- list()

  for(i in 1:nRecords){
    temp.curve <- records[[i]]

    temp.recordType <- temp.curve@recordType
    temp.curveType <- temp.curve@curveType

    temp.metadata <- temp.curve@info

    OLD_DTYPE<- temp.metadata$DTYPE
    temp.metadata<-cbind(temp.metadata,OLD_DTYPE)

    temp.temperatures <- as.numeric(temp.curve@data[,1])
    temp.data <- as.numeric(temp.curve@data[,2])
    temp.error <- abs(temp.curve*relative.error)

    temp.analysis <- list()
    temp.reserved <- list()

    temp.TLum.data.curve <- set_TLum.Data.Curve(recordType = temp.recordType,
                                                curveType = temp.curveType,
                                                temperatures= temp.temperatures,
                                                data = temp.data,
                                                error = temp.error,
                                                metadata = temp.metadata,
                                                analysis = temp.analysis,
                                                .RESERVED = temp.reserved)

    new.records <- c(new.records, temp.TLum.data.curve)
  }

  new.TLum.Analysis <- set_TLum.Analysis(new.records, new.protocol)

  return(new.TLum.Analysis)
}
