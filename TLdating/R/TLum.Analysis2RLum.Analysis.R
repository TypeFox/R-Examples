#' Convert RLum.Analysis object to an TLum.Analysis.
#'
#' This function convert a \linkS4class{TLum.Analysis} object into a \linkS4class{RLum.Analysis} object, from the 'Luminescence' package.
#' The 'TLdating' package is implemented to only works with its own class of object (TLum.Analysis, TLum.Analysis and TLum.Data.Curve).
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the TL curves used for the ED estimation.
#'
#' @details
#' This function use the data from the TLum.Analysis to create a new RLum.Analysis.
#' During the process, all information relative to the uncertainties and stored in the TLum.Analysis object are lost.
#' The original data-type of each luminescence curve is also restored.
#'
#'
#' @author David Strebler, University of Cologne (Germany), \cr David Strebler
#'
#' @export TLum.Analysis2RLum.Analysis


TLum.Analysis2RLum.Analysis <- function(
  object

){

  # Integrity Check ---------------------------------------------------------
  if(missing(object)){
    stop("[TLum.Analysis2RLum.Analysis] Error: Input object is missing.")

  }else if(!is(object,"TLum.Analysis")){
    stop("[TLum.Analysis2RLum.Analysis] Error: Input object is not of type 'RLum.Analysis'.")
  }

  # ------------------------------------------------------------------------------

  nRecords <- length(object@records)

  protocol <- object@protocol
  records <- object@records

  new.records <- list()

  for(i in 1:nRecords){
    temp.curve <- records[[i]]

    temp.recordType <- temp.curve@recordType
    temp.curveType <- temp.curve@curveType

    temp.data <- cbind(temp.curve@temperatures,temp.curve@data)

    temp.info <- temp.curve@metadata

    temp.info$DTYPE <- temp.info$OLD_DTYPE
    temp.info$OLD_DTYPE <- NULL

    temp.originator <- "TLum.Analysis2RLum.Analysis"

    # temp.RLum.Data.Curve <- set_RLum.Data.Curve(recordType=temp.recordType,
    #                                             curveType=temp.curveType,
    #                                             data=temp.data,
    #                                             info=temp.info)

    temp.RLum.Data.Curve <- new("RLum.Data.Curve")
    temp.RLum.Data.Curve@recordType <- temp.recordType
    temp.RLum.Data.Curve@curveType <- temp.curveType
    temp.RLum.Data.Curve@data <- temp.data
    temp.RLum.Data.Curve@info <- temp.info
    temp.RLum.Data.Curve@originator <- temp.originator

    new.records <- c(new.records, temp.RLum.Data.Curve)
  }

  # new.RLum.Analysis <- set_RLum.Analysis(records=new.records,
  #                                        protocol=protocol)

  new.RLum.Analysis <- new("RLum.Analysis")
  new.RLum.Analysis@records <- new.records
  new.RLum.Analysis@protocol <- protocol

  return(new.RLum.Analysis)
}
