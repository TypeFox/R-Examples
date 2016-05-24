#' Convert Risoe.BINfileData object to an TLum.BIN.File object.
#'
#' This function convert a \linkS4class{Risoe.BINfileData} object, created using the package'Luminescence', into a TLum.BIJ.file object, which will be used by this package.
#' The TL package is implemented to only works with its own class of object (TLum.BIN.File, TLum.Analysis and TLum.Data.Curve).
#'
#' @param object
#'  \code{\linkS4class{Risoe.BINfileData}} (\bold{required}): object containing the TL curves used for the ED estimation.
#' @param relative.error
#'  \link{numeric} (\bold{required}): Relative error of the TL signals. Generally, it is between 0.02 and 0.1.
#'
#' @details
#' This function use the data from the Risoe.BINFileData and the relative.error specified to create a absolute error matrix.
#' Then it create a new TLum.BIN.File including all the information from the Risoe.BINFileData and the new absolute error matrix.
#' For practical reason, the TLdating package considers the error as random. It means that the systematic component of the error will be ignored.
#'
#' @seealso \link{TLum.BIN.File2Risoe.BINfileData}, \link{TLum.BIN.File2TLum.Analysis} and \link{TLum.BIN.File2TLum.Data.Curve}.
#'
#' @author David Strebler, University of Cologne (Germany), \cr David Strebler
#'
#' @export Risoe.BINfileData2TLum.BIN.File


Risoe.BINfileData2TLum.BIN.File <- function(
  object,

  relative.error

){

  # Integrity Check ---------------------------------------------------------
  if(missing(object)){
    stop("[Risoe.BINfileData2TLum.BIN.File] Error: Input object is missing.")

  }else if(!is(object,"Risoe.BINfileData")){
    stop("[Risoe.BINfileData2TLum.BIN.File] Error: Input object is not of type 'Risoe.BINfileData'.")
  }

  if(missing(relative.error)){
    stop("[Risoe.BINfileData2TLum.BIN.File] Error: Input relative error is missing.")

  }else if(!is.numeric(relative.error)){
    stop("[Risoe.BINfileData2TLum.BIN.File] Error: Relative error is not of type 'numeric'.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Value check
  if(relative.error > 1){
    warning("[Risoe.BINfileData2TLum.BIN.File] Warning: Input 'relative.error' > 1.")

  }else if(relative.error < -1){
    relative.error <- abs(relative.error)
    warning("[Risoe.BINfileData2TLum.BIN.File] Warning: Input 'relative.error' < -1.")

  }else if(relative.error < 0){
    relative.error <- abs(relative.error)
    warning("[Risoe.BINfileData2TLum.BIN.File] Warning: Input 'relative.error' < 0.")
  }
  # ------------------------------------------------------------------------------

  # METADATA
  metadata <- object@METADATA
  OLD_DTYPE<- metadata$DTYPE

  new.metadata<-cbind(metadata,OLD_DTYPE)

  # DATA
  new.data <- object@DATA

  # ERROR
  nRecord <- length(new.data)

  new.error <- list()

  for(i in 1:nRecord){

    temp.curve <- new.data[[i]]
    temp.absolute.error <- abs(temp.curve*relative.error)
    temp.error <- list(temp.absolute.error)

    new.error <- c(new.error, temp.error)

  }

  # .RESERVED
  new.reserved <- object@.RESERVED

  # FileData
  new.TLum.BIN.File <- set_TLum.BIN.File(METADATA = new.metadata,
                                         DATA = new.data,
                                         ERROR = new.error,
                                         .RESERVED = new.reserved)

  return(new.TLum.BIN.File)
}
