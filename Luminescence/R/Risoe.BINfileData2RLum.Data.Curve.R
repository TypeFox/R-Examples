#' Convert an element from a Risoe.BINfileData object to an RLum.Data.Curve
#' object
#'
#' The function converts one specified single record from a Risoe.BINfileData
#' object to an RLum.Data.Curve object.
#'
#' The function extracts all \code{METADATA} from the \code{Risoe.BINfileData}
#' object and stores them in the \code{RLum.Data.Curve} object. This function
#' can be used stand-alone, but is the base function for \code{\link{Risoe.BINfileData2RLum.Analysis}}.
#'
#' @param object \code{\linkS4class{Risoe.BINfileData}} (\bold{required}):
#' \code{Risoe.BINfileData} object
#'
#' @param id \code{\link{integer}} (\bold{required}): record id in the
#' \code{Risoe.BINfileData} object of the curve that is to be stored in the
#' \code{RLum.Data.Curve} object. If no value for id is provided, the record
#' has to be specified by \code{pos}, \code{set} and \code{run}.
#'
#' @param pos \code{\link{integer}} (optional): record position number in the
#' \code{Risoe.BINfileData} object of the curve that is to be stored in the
#' \code{RLum.Data.Curve} object. If a value for \code{id} is provided, this
#' argument is ignored.
#'
#' @param run \code{\link{integer}} (optional): record run number in the
#' \code{Risoe.BINfileData} object of the curve that is to be stored in the
#' \code{RLum.Data.Curve} object. If a value for \code{id} is provided, this
#' argument is ignored.
#'
#' @param set \code{\link{integer}} (optional): record set number in the
#' \code{Risoe.BINfileData} object of the curve that is to be stored in the
#' \code{RLum.Data.Curve} object. If a value for \code{id} is provided, this
#' argument is ignored.
#'
#' @return Returns an \code{\linkS4class{RLum.Data.Curve}} object.
#'
#' @note Due to changes in the BIN-file (version 3 to version 4) format the recalculation of TL-curves might be not
#' overall correct for cases where the TL measurement is combined with a preheat.
#'
#' @section Function version: 0.2.1
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France),
#' Christoph Burow, Universtiy of Cologne (Germany)
#'
#' @seealso \code{\link{Risoe.BINfileData2RLum.Analysis}},
#' \code{\link{set_RLum}}, \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Analysis}}, \code{\linkS4class{Risoe.BINfileData}},
#' \code{\link{plot_RLum}}
#'
#' @references #
#'
#' @keywords manip
#'
#' @examples
#'
#'
#' ##get package example data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##convert one record
#' Risoe.BINfileData2RLum.Data.Curve(CWOSL.SAR.Data, id = 1)
#'
#' @export
Risoe.BINfileData2RLum.Data.Curve <- function(
  object,
  id,
  pos,
  run,
  set
){


  # Integrity Check ---------------------------------------------------------

  if (!is(object,"Risoe.BINfileData")){
    stop("[Risoe.BINfileData2RLum.Data.Curve()] Input object is not of type 'Risoe.BINfileData'.")
  }

  ##if id is set, no input for pos and rund is nescessary
  if(missing(id)){

    if(missing(pos) == TRUE | missing(run) == TRUE | missing(set) == TRUE){

      temp.missing.arguments <- paste(c(if(missing(pos)==TRUE){"pos"},
                                        if(missing(set)==TRUE){"set"},
                                        if(missing(run)==TRUE){"run"}), collapse=", ")

      stop(paste("[Risoe.BINfileData2RLum.Data.Curve()] Arguments are missing: ",
                 temp.missing.arguments, ". Or set id.", sep = ""))

    }

    if (is(pos,"numeric")==FALSE){
      stop("[Risoe.BINfileData2RLum.Data.Curve()] Argument 'pos' has to be of data type integer.")
    }

    if (is(set,"numeric")==FALSE){
      stop("[Risoe.BINfileData2RLum.Data.Curve()]Argument 'set' has to be of data type integer.")
    }

    if (is(run,"numeric")==FALSE){
      stop("[Risoe.BINfileData2RLum.Data.Curve()] Argument 'run' has to be of data type integer.")
    }

    if (length(which(pos/1:48 == 1)) == 0){
      stop("[Risoe.BINfileData2RLum.Data.Curve()] Value for 'pos' out of bounds.")
    }

    ##get and check valid positions
    positions.valid <- paste(as.character(unique(object@METADATA[,"POSITION"])), collapse=", ")

    if ((pos %in% unique(object@METADATA[,"POSITION"])) == FALSE){
      stop(paste("[Risoe.BINfileData2RLum.Data.Curve] Error: pos = ",pos, " is not valid.
               Valid positions are: ", positions.valid, sep=""))
    }

    ##get and check valid positions
    positions.valid <- paste(as.character(unique(object@METADATA[,"SET"])), collapse=", ")

    if ((set %in% unique(object@METADATA[,"SET"])) == FALSE){
      stop(paste("[Risoe.BINfileData2RLum.Data.Curve] Error: set = ",set, " is not valid.
               Valid values are: ", positions.valid, sep=""))
    }


    ##get and check valid positions
    positions.valid <- paste(as.character(unique(object@METADATA[,"RUN"])), collapse=", ")

    if ((run %in% unique(object@METADATA[,"RUN"])) == FALSE){
      stop(paste("[Risoe.BINfileData2RLum.Data.Curve] Error: run = ",run, " is not valid.
               Valid values are: ", positions.valid, sep=""))
    }

  }else{

    ##check if id is valid at all
    temp.range.id <- range(object@METADATA[,"ID"])

    if ((id %in% unique(object@METADATA[,"ID"])) == FALSE){
      stop(paste("[Risoe.BINfileData2RLum.Data.Curve()] id = ",id, " is not a valid record id. Allowed value range ", min(temp.range.id), " : ", max(temp.range.id),".", sep=""))

    }

  }


  # grep id of record -------------------------------------------------------

  ##if id is set, no input for pos and rund is nescessary
  if (missing(id)) {
    id <- object@METADATA[object@METADATA[, "POSITION"] == pos &
                            object@METADATA[, "SET"] == set &
                            object@METADATA[, "RUN"] == run,
                          "ID"]

  }


  # Select values -----------------------------------------------------------

  ##build matrix
  if(object@METADATA[id,"NPOINTS"][1] != 0){

    ##set variables
    temp.x <- vector(mode = "numeric", length = object@METADATA[id,"NPOINTS"])
    temp.y <- vector(mode = "integer", length = object@METADATA[id,"NPOINTS"])

    if(object@METADATA[id, "LTYPE"] == "TL" && as.numeric(object@METADATA[id, "VERSION"]) >=4){

      temp.x <- c(
        seq(
          from = object@METADATA[id, "LOW"],
          to = object@METADATA[id, "AN_TEMP"],
          length.out = object@METADATA[id, "TOLDELAY"]
        ),
        seq(
          from = object@METADATA[id, "AN_TEMP"],
          to = object@METADATA[id, "AN_TEMP"],
          length.out = object@METADATA[id, "TOLON"]
        ),
        seq(
          from = object@METADATA[id, "AN_TEMP"],
          to = object@METADATA[id, "HIGH"],
          length.out = object@METADATA[id, "TOLOFF"]
        )
      )

    }else{

      temp.x <- seq(
        from = object@METADATA[id, "LOW"],
        to = object@METADATA[id, "HIGH"],
        length.out = object@METADATA[id, "NPOINTS"]
      )

    }

    temp.y <- unlist(object@DATA[id])


  }else{
    temp.x <- NA
    temp.y <- NA

    warning("NPOINTS was 0, RLum.Data.Curve-object with NA-values produced.")

  }


  temp.data <- matrix(c(temp.x,temp.y), ncol=2, byrow=FALSE)
  temp.recordType <- as.character(object@METADATA[id,"LTYPE"])
  temp.info <- as.list(object@METADATA[id,])


  # Build object ------------------------------------------------------------

  newRLumDataCurve.Risoe.BINfileData2RLum.Data.Curve <- set_RLum(
    class = "RLum.Data.Curve",
    recordType = temp.recordType,
    data = temp.data,
    info = temp.info)

  return(newRLumDataCurve.Risoe.BINfileData2RLum.Data.Curve)

}
