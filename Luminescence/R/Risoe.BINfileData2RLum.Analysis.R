#' Convert Risoe.BINfileData object to an RLum.Analysis object
#'
#' Converts values from one specific position of a Risoe.BINfileData S4-class
#' object to an RLum.Analysis object. Internally the function uses the function
#' \code{\link{Risoe.BINfileData2RLum.Data.Curve}} to recalculate the curves.
#'
#' The \code{\linkS4class{RLum.Analysis}} object requires a set of curves for
#' specific further protocol analyses. However, the
#' \code{\linkS4class{Risoe.BINfileData}} usually contains a set of curves for
#' different aliquots and different protocol types that may be mixed up.
#' Therefore, a conversion is needed.
#'
#' @param object \code{\linkS4class{Risoe.BINfileData}} (\bold{required}):
#' \code{Risoe.BINfileData} object
#'
#' @param pos \code{\link{numeric}} (optional): position number of the
#' \code{Risoe.BINfileData} object for which the curves are stored in the
#' \code{RLum.Analysis} object. If \code{length(position)>1} a list of \code{RLum.Analysis} objects
#' is returned. If nothing is provided every position will be converted. If the position is not valid \code{NA} is
#' returned.
#'
#' @param grain \code{\link{vector}, \link{numeric}} (optional): grain number from
#' the measurement to limit the converted data set (e.g., \code{grain =
#' c(1:48)}). Please be aware that this option may lead to unwanted effects, as the output
#' is strictly limited to the choosen grain number for all position numbers
#'
#' @param run \code{\link{vector}, \link{numeric}} (optional): run number from
#' the measurement to limit the converted data set (e.g., \code{run =
#' c(1:48)}).
#'
#' @param set \code{\link{vector}, \link{numeric}} (optional): set number from
#' the measurement to limit the converted data set (e.g., \code{set =
#' c(1:48)}).
#'
#' @param ltype \code{\link{vector}, \link{character}} (optional): curve type
#' to limit the converted data. Allowed values are: \code{IRSL}, \code{OSL},
#' \code{TL}, \code{RIR}, \code{RBR} and \code{USER}
#'
#' @param protocol \code{\link{character}} (optional): sets protocol type for
#' analysis object. Value may be used by subsequent analysis functions.
#'
#' @param txtProgressBar \link{logical} (with default): enables or disables
#' \code{\link{txtProgressBar}}.
#'
#' @return Returns an \code{\linkS4class{RLum.Analysis}} object.
#'
#' @note The \code{protocol} argument of the \code{\linkS4class{RLum.Analysis}}
#' object is set to 'unknown' if not stated otherwise.
#'
#' @section Function version: 0.2.2
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @seealso \code{\linkS4class{Risoe.BINfileData}}, \code{\link{Risoe.BINfileData2RLum.Data.Curve}}
#' \code{\linkS4class{RLum.Analysis}}, \code{\link{read_BIN2R}}
#'
#' @references #
#'
#' @keywords manip
#'
#' @examples
#'
#' ##load data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##convert values for position 1
#' Risoe.BINfileData2RLum.Analysis(CWOSL.SAR.Data, pos = 1)
#'
#' @export
Risoe.BINfileData2RLum.Analysis<- function(
  object,
  pos,
  grain,
  run,
  set,
  ltype,
  protocol = "unknown",
  txtProgressBar = FALSE
){


  # Integrity Check ---------------------------------------------------------

  if (!is(object,"Risoe.BINfileData")){
    stop("[Risoe.BINfileData2RLum.Analysis()] Input object is not of type 'Risoe.BINfileData'.")
  }

  if (missing(pos)){
    pos <- unique(object@METADATA[,"POSITION"])
  }

  if (!is(pos,"numeric")){
    stop("[Risoe.BINfileData2RLum.Analysis] Argument 'pos' has to be of type numeric.")
  }

  ##get and check valid positions
  positions.valid <- paste(as.character(unique(object@METADATA[,"POSITION"])), collapse=", ")


  if (!all(pos %in% unique(object@METADATA[,"POSITION"]))){

    warning(paste("[Risoe.BINfileData2RLum.Analysis] Error: pos=",pos, " invalid.
              Valid positions are: ", positions.valid, sep=""))

    ##flag position
    pos.valid <- FALSE

  }else{

    pos.valid <- TRUE

  }

  ##WARNINGS
  if (length(which(max(pos)/1:48 == 1)) == 0){
    warning("[Risoe.BINfileData2RLum.Analysis] Value for 'pos' out bounds specified for
            a Risoe BIN-file.")
  }


  # Grep run and set data ---------------------------------------------------

  if(pos.valid){
    ##grep values according to their criteria and check for validity

    ##grain
    if (missing(grain)) {
      grain <- unique(object@METADATA[, "GRAIN"])
    }
    else{

      if(grain %in% unique(unique(object@METADATA[, "GRAIN"])) == FALSE){

        ##get only valid grain numbers
        grain.valid <- paste(as.character(unique(object@METADATA[,"GRAIN"])), collapse=", ")

        stop(paste("[Risoe.BINfileData2RLum.Analysis()] grain = ", grain, " contain invalid run(s).
                   Valid grain values are: ", grain.valid, sep=""))


      }


    }

    ##run
    if(missing(run)){run <- unique(object@METADATA[, "RUN"])} else{

      if(TRUE %in% unique(unique(object@METADATA[, "RUN"]) %in% run) != TRUE){

        ##get and check valid positions
        run.valid <- paste(as.character(unique(object@METADATA[,"RUN"])), collapse=", ")

        stop(paste("[Risoe.BINfileData2RLum.Analysis()] run = ", run, " contain invalid run(s).
                   Valid runs are: ", run.valid, sep=""))

      }

    }

    #set
    if(missing(set)){set <- unique(object@METADATA[, "SET"])} else{

      if(TRUE %in% unique(unique(object@METADATA[, "SET"]) %in% set) != TRUE){

        ##get and check valid positions
        set.valid <- paste(as.character(unique(object@METADATA[,"SET"])), collapse=", ")

        stop(paste("[Risoe.BINfileData2RLum.Analysis] set = ", set, " contain invalid set(s).
                   Valid sets are: ", set.valid, sep=""))

      }

    }

    ##ltype
    if(missing(ltype)){ltype <- unique(object@METADATA[, "LTYPE"])} else{

      if(TRUE %in% unique(unique(object@METADATA[, "LTYPE"]) %in% ltype) != TRUE){

        ##get and check valid positions
        ltype.valid <- paste(as.character(unique(object@METADATA[,"LTYPE"])), collapse=", ")

        stop(paste("[Risoe.BINfileData2RLum.Analysis] ltype = ", ltype, " contain invalid ltype(s).
               Valid ltypes are: ", ltype.valid, sep=""))

      }

    }


    # Select values and convert them-----------------------------------------------------------

    ##set progressbar to false if only one position is provided
    if(txtProgressBar & length(pos)<2){
      txtProgressBar <- FALSE

    }

      ##set progress bar
      if(txtProgressBar){
        pb<-txtProgressBar(min=min(pos),max=max(pos), char="=", style=3)
      }

    object <- lapply(pos, function(pos){

      ##update progress bar
      if(txtProgressBar==TRUE){
        setTxtProgressBar(pb, value = pos)
      }

      ##deselect all values
      object@METADATA[, "SEL"] <- FALSE

      ##select data
      object@METADATA[which(
        object@METADATA[,"POSITION"] == pos &
          object@METADATA[,"GRAIN"] %in% grain &
          object@METADATA[,"RUN"] %in% run &
          object@METADATA[,"SET"] %in% set &
          object@METADATA[,"LTYPE"] %in% ltype
      )
      , "SEL"] <- TRUE

      # Limit object to selection -----------------------------------------------

      object@DATA <-
        object@DATA[object@METADATA[object@METADATA[,"SEL"],"ID"]]

      object@METADATA <-
        object@METADATA[object@METADATA[,"SEL"],]

      ##correct ID values after limitation, if we don't do that we get problems with
      ##the conversion later on
      object@METADATA$ID <- 1:nrow(object@METADATA)


      # Convert values ----------------------------------------------------------
      object <- set_RLum(
        class = "RLum.Analysis",
        records = lapply(object@METADATA$ID,function(x) {
          Risoe.BINfileData2RLum.Data.Curve(object, id = x)
        }),
        protocol = protocol,
        originator = "Risoe.BINfileData2RLum.Analysis"
      )

      return(object)
    })

    if(txtProgressBar){close(pb)}

    ##this is necassary to not break with previous code
    if(length(object) == 1){
      invisible(object[[1]])

    }else{
      invisible(object)

    }


  }else{

    invisible(NA)

  }
}
