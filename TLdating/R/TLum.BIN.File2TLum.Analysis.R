#' Convert TLum.BIN.File object to an TLum.Analysis object.
#'
#' This function convert a \linkS4class{TLum.BIN.File} in a \linkS4class{TLum.Analysis} object.
#' A \linkS4class{TLum.Analysis} object is a list of \linkS4class{TLum.Data.Curve} object.
#' It is possible to specify which luminescence curves will be keeped.
#'
#' @param object
#'  \code{\linkS4class{TLum.BIN.File}} (\bold{required}): object containing the luminescence curves.
#' @param protocol
#'  \link{character} (with default): protocol used.
#' @param rec_ramp2PH
#'  \link{logical} (with default): Indicate if the signal was record during the ramp up to the preheat temperature.
#' @param rec_duringPH
#'  \link{logical} (with default): Indicate if the signal was record during the preheat plateau.
#'
#' @return
#'  This function will return a \linkS4class{TLum.Analysis} object.
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


TLum.BIN.File2TLum.Analysis <- function(

  object,
  protocol = "unknown",
  rec_duringPH =TRUE,
  rec_ramp2PH =TRUE

){

  # Integrity Check ---------------------------------------------------------

  if (!is(object,"TLum.BIN.File")){
    stop("[TLum.BIN.File] Error: Input object is not of type 'TLum.BIN.File'.")
  }

  if(!is.logical(rec_ramp2PH) || is.na(rec_ramp2PH)){
    stop("[calc_TL.temperature] Error: Input 'rec_ramp2PH' is not of type 'logical'.")
  }

  if(!is.logical(rec_duringPH) || is.na(rec_duringPH)){
    stop("[calc_TL.temperature] Error: Input 'rec_duringPH' is not of type 'logical'.")
  }

  if(!is.character(protocol)){
    stop("[calc_TL.temperature] Error: Input 'protocol' is not of type 'character'.")
  }

  # ------------------------------------------------------------------------------

  nRecords <- length(object@DATA)

  positions <- as.integer(object@METADATA$POSITION)
  runs <- as.integer(object@METADATA$RUN)
  sets <- as.integer(object@METADATA$SET)
  ltypes <- as.character(object@METADATA$LTYPE)
  ids <- as.integer(object@METADATA$ID)

  # ------------------------------------------------------------------------------
  # value check


  new.records <-list ()

  for(i in 1 : nRecords){

    temp.pos <- positions[i]
    temp.run <- runs[i]
    temp.set <- sets[i]
    temp.ltype <- ltypes[i]
    temp.id <- ids[i]

    temp.record <- TLum.BIN.File2TLum.Data.Curve(object = object,
                                                 pos = temp.pos,
                                                 run = temp.run,
                                                 set = temp.set,
                                                 id = temp.id,
                                                 rec_duringPH =rec_duringPH,
                                                 rec_ramp2PH =rec_ramp2PH)

    new.records <- c(new.records, temp.record)
  }

  new.TLum.Analysis <- set_TLum.Analysis(records = new.records,
                                         protocol = protocol)

  return(new.TLum.Analysis)
}
