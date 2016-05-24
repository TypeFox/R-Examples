#' Convert TLum.Analysis object to an TLum.BIN.File object.
#'
#' This function convert a \linkS4class{TLum.BIN.File} in a \linkS4class{TLum.Analysis} object.
#' A \linkS4class{TLum.Analysis} object is a list of \linkS4class{TLum.Data.Curve} object.
#' It is possible to specify which luminescence curves will be keeped.
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the luminescence curves.
#'
#' @return
#'  This function will return a \linkS4class{TLum.BIN.File} object.
#'
#' @seealso
#'  \linkS4class{TLum.Analysis},
#'  \linkS4class{TLum.BIN.File},
#'  \linkS4class{TLum.Data.Curve} and
#'  \link{TLum.BIN.File2TLum.Data.Curve}.
#'
#' @author David Strebler, University of Cologne (Germany), \cr David Strebler
#'
#' @export TLum.BIN.File2TLum.Analysis


TLum.Analysis2TLum.BIN.File <- function(

  object
){

  # Integrity Check ---------------------------------------------------------

  if (!is(object,"TLum.Analysis")){
    stop("[TLum.Analysis2TLum.BIN.File] Error: Input object is not of type 'TLum.Analysis'.")
  }

  # ------------------------------------------------------------------------------

  nRecords <- length(object@records)

  metadata.names <- names(object@records[[1]]@metadata)

  new.METADATA <- data.frame(row.names = metadata.names)
  new.DATA <- list()
  new.ERROR <- list()
  new.RESERVED <- list()


  for(i in 1 : nRecords){

    temp.curves <- object@records[[i]]

    temp.metadata <- as.data.frame(get_TLum.Data.Curve(object = temp.curves,ref = "metadata"))
    temp.data <- list(get_TLum.Data.Curve(object = temp.curves,ref = "data"))
    temp.error <- list(get_TLum.Data.Curve(object = temp.curves,ref = "error"))
    temp.reserved <- list(get_TLum.Data.Curve(object = temp.curves,ref = ".RESERVED"))


    new.METADATA <- rbind(new.METADATA, temp.metadata)
    new.DATA <- c(new.DATA, temp.data)
    new.ERROR <- c(new.ERROR, temp.error)
    new.RESERVED <- c(new.RESERVED, temp.reserved)



  }

  new.TLum.BIN.File <- set_TLum.BIN.File(METADATA = new.METADATA,
                                         DATA = new.DATA,
                                         ERROR = new.ERROR,
                                         .RESERVED = new.RESERVED)

  return(new.TLum.BIN.File)
}
