#' Convert an element from a \linkS4class{TLum.BIN.File} object into a \linkS4class{TLum.Data.Curve} objet
#'
#' This function extract a curve from a \linkS4class{TLum.BIN.File} object and convert it into a \linkS4class{TLum.Data.Curve} objet.
#' The extract element can be identify either by its id or by its position, run and set.
#'
#' @param object
#'  \code{\linkS4class{TLum.BIN.File}} (\bold{required}): object containing the luminescence curves.
#' @param id
#'  \link{integer} (with default): id of the curve.
#' @param pos
#'  \link{integer} (with default): position of the curve.
#' @param run
#'  \link{integer} (with default): run of the curve.
#' @param set
#'  \link{integer} (with default): set of the curve.
#' @param rec_ramp2PH
#'  \link{logical} (with default): Indicate if the signal was record during the ramp up to the preheat temperature.
#' @param rec_duringPH
#'  \link{logical} (with default): Indicate if the signal was record during the preheat plateau.
#'
#' @details
#'  The element that is extracted to be converted into a \linkS4class{TLum.Data.Curve} objet can be identify
#'  either by its id or by its position, run and set.
#'
#' @return
#'  This function return a \linkS4class{TLum.Data.Curve} objet.
#'
#' @export TLum.BIN.File2TLum.Data.Curve
#'

TLum.BIN.File2TLum.Data.Curve <- function(
  object,
  id,
  pos,
  run,
  set,
  rec_duringPH =TRUE,
  rec_ramp2PH =TRUE

){


  # Integrity Check ---------------------------------------------------------

  if (is(object,"TLum.BIN.File")==FALSE){
    stop("[TLum.BIN.File2TLum.Data.Curve] Error: Input object is not of type 'TLum.BIN.File'.")
  }

  if(!is.logical(rec_ramp2PH) || is.na(rec_ramp2PH)){
    stop("[calc_TL.temperature] Error: Input 'rec_ramp2PH' is not of type 'logical'.")
  }

  if(!is.logical(rec_duringPH) || is.na(rec_duringPH)){
    stop("[calc_TL.temperature] Error: Input 'rec_duringPH' is not of type 'logical'.")
  }

  ##if id is set, no input for pos and run is nescessary
  if(missing(id) == TRUE){

    if(missing(pos) == TRUE | missing(run) == TRUE | missing(set) == TRUE){

      temp.missing.arguments <- paste(c(if(missing(pos)==TRUE){"pos"},
                                        if(missing(set)==TRUE){"set"},
                                        if(missing(run)==TRUE){"run"}), collapse=", ")

      stop(paste("[TLum.BIN.File2TLum.Data.Curve] Error: Arguments are missing: ",
                 temp.missing.arguments, ". Or set id.", sep = ""))

    }

    if (is(pos,"numeric")==FALSE){
      stop("[TLum.BIN.File2TLum.Data.Curve] Error: Argument 'pos' has to be of data type integer.")
    }

    if (is(set,"numeric")==FALSE){
      stop("[TLum.BIN.File2TLum.Data.Curve] Error: Argument 'set' has to be of data type integer.")
    }

    if (is(run,"numeric")==FALSE){
      stop("[TLum.BIN.File2TLum.Data.Curve] Error: Argument 'run' has to be of data type integer.")
    }

    if (length(which(pos/1:48 == 1)) == 0){
      stop("[TLum.BIN.File2TLum.Data.Curve] Error: Value for 'pos' out of bounds.")
    }

    ##get and check valid positions
    positions.valid <- paste(as.character(unique(object@METADATA[,"POSITION"])), collapse=", ")

    if ((pos %in% unique(object@METADATA[,"POSITION"])) == FALSE){
      stop(paste("[TLum.BIN.File2TLum.Data.Curve] Error: pos = ",pos, " is not valid.
                 Valid positions are: ", positions.valid, sep=""))
    }

    ##get and check valid positions
    positions.valid <- paste(as.character(unique(object@METADATA[,"SET"])), collapse=", ")

    if ((set %in% unique(object@METADATA[,"SET"])) == FALSE){
      stop(paste("[TLum.BIN.File2TLum.Data.Curve] Error: set = ",set, " is not valid.
                 Valid values are: ", positions.valid, sep=""))
    }


    ##get and check valid positions
    positions.valid <- paste(as.character(unique(object@METADATA[,"RUN"])), collapse=", ")

    if ((run %in% unique(object@METADATA[,"RUN"])) == FALSE){
      stop(paste("[TLum.BIN.File2TLum.Data.Curve] Error: run = ",run, " is not valid.
                 Valid values are: ", positions.valid, sep=""))
    }

    }else{

      ##check if id is valid
      temp.range.id <- range(object@METADATA[,"ID"])

      if ((id %in% unique(object@METADATA[,"ID"])) == FALSE){
        stop(paste("[TLum.BIN.File2TLum.Data.Curve] Error: id = ",id, " is not a valid record id. Allowed value range ", min(temp.range.id), " : ", max(temp.range.id),".", sep=""))

      }

    }


  # grep id of record -------------------------------------------------------

  ##if id is set, no input for pos and rund is nescessary
  if(missing(id) == TRUE){

    id <- object@METADATA[object@METADATA[,"POSITION"] == pos &
                            object@METADATA[,"SET"] == set &
                            object@METADATA[,"RUN"] == run,
                          "ID"]


  }


  # Select values -----------------------------------------------------------

  ##build matrix
  Tmax <- object@METADATA[id,"HIGH"]
  nPoints <- object@METADATA[id,"NPOINTS"]
  Hrate <- object@METADATA[id,"RATE"]
  an_time  <- object@METADATA[id,"AN_TIME"]
  an_temp  <- object@METADATA[id,"AN_TEMP"]

  temperatures.data <- calc_TL.temperature(nPoints = nPoints,
                                          Tmax = Tmax,
                                          Hrate = Hrate,
                                          an_temp = an_temp,
                                          an_time = an_time,
                                          rec_ramp2PH = rec_ramp2PH,
                                          rec_duringPH = rec_duringPH)

  new.temperatures <- get_TLum.Results(temperatures.data,"temperatures")

  new.data <- as.numeric(unlist(object@DATA[id]))

  new.error <- as.numeric(unlist(object@ERROR[id]))

  new.recordType <- as.character(object@METADATA[id,"LTYPE"])

  new.metadata <- as.list(object@METADATA[id,])

  new.analysis <- list()

  new.RESERVED <- object@.RESERVED[[id]]

  # Build object ------------------------------------------------------------

  new.TLum.Data.Curve <- set_TLum.Data.Curve(recordType = new.recordType,
                                             temperatures = new.temperatures,
                                             data = new.data,
                                             error = new.error,
                                             metadata = new.metadata,
                                             analysis = new.analysis,
                                             .RESERVED = new.RESERVED)

  return(new.TLum.Data.Curve)
    }
