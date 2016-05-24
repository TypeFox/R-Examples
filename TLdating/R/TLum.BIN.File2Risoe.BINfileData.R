#' Convert TLum.BIN.File object to an Risoe.BINfileData object.
#'
#' This function convert \linkS4class{TLum.BIN.File} object into a \linkS4class{Risoe.BINfileData} object that is usable by the \link{Luminescence} package.
#'
#' @param object
#'  \code{\linkS4class{TLum.BIN.File}} (\bold{required}): object containing the TL curves used for the ED estimation.
#'
#' @return
#'  This function return an \linkS4class{Risoe.BINfileData} containing all information previously stored in the \linkS4class{TLum.BIN.File} except the uncertainties matrix.
#'  To avoid conflicts with other software, the original data type of each curves is restored.
#'
#' @author David Strebler
#'
#' @export TLum.BIN.File2Risoe.BINfileData

TLum.BIN.File2Risoe.BINfileData <- function(

  object
){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(object)){
    stop("[TLum.BIN.File2Risoe.BINfileData] Error: Input 'object' is missing.")
  }else if (!is(object,"TLum.BIN.File")){
    stop("[TLum.BIN.File2Risoe.BINfileData] Error: Input 'object' is not of type 'TLum.BIN.File'.")
  }
  # ------------------------------------------------------------------------------

  # METADATA
  metadata <- object@METADATA
  old.dtype <- metadata$OLD_DTYPE

  new.metadata <- metadata[names(metadata)!="OLD_DTYPE"]
  new.metadata$DTYPE <- old.dtype

  new.metadata$COMMENT <- as.character(new.metadata$COMMENT)

  # DATA
  new.data <- object@DATA

  # .RESERVED
  new.reserved <- object@.RESERVED

  # Risoe.BINfileData
  new.bin <- new("Risoe.BINfileData",
                 METADATA=new.metadata,
                 DATA=new.data,
                 .RESERVED = new.reserved)

  return(new.bin)
}
