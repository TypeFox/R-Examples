#' Extract irradiation times from an XSYG file
#'
#' Extracts irradiation times, dose and times since last irradiation, from a
#' Freiberg Instruments XSYG-file. These information can be further used to
#' update an existing BINX-file
#'
#' The function was written to compensate missing information in the BINX-file
#' output of Freiberg Instruments lexsyg readers. As all information are
#' available within the XSYG-file anyway, these information can be extracted
#' and used for further analysis or/and to stored in a new BINX-file, which can
#' be further used by other software, e.g. Analyst (Geoff Duller). \cr
#'
#' Typical application example: g-value estimation from fading measurements
#' using the Analyst or any other self written script.\cr
#'
#' Beside the some simple data transformation steps the function applies the
#' functions \code{\link{read_XSYG2R}}, \code{\link{read_BIN2R}},
#' \code{\link{write_R2BIN}} for data import and export.
#'
#' @param object \code{\link{character}} (\bold{required}) or
#' \code{\linkS4class{RLum.Analysis}} object: path and file name of the XSYG
#' file or an \code{\linkS4class{RLum.Analysis}} produced by the function
#' \code{\link{read_XSYG2R}}. \cr
#'
#' \bold{Note}: If an \code{\linkS4class{RLum.Analysis}} is used, any input for
#' the arguments \code{file.BINX} and \code{recordType} will be ignored!
#' @param file.BINX \code{\link{character}} (optional): path and file name of
#' an existing BINX-file. If a file name is provided the file will be updated
#' with the information from the XSYG file in the same folder as the original
#' BINX-file.\cr Note: The XSYG and the BINX-file have to be originate from the
#' same measurement!
#' @param recordType \code{\link{character}} (with default): select relevant
#' curves types from the XSYG file or \code{\linkS4class{RLum.Analysis}}
#' object. As the XSYG-file format comprises much more information than usually
#' needed for routine data analysis and allowed in the BINX-file format, only
#' the relevant curves are selected by using the function
#' \code{\link{get_RLum}}. The argument \code{recordType} works as
#' described for this function. \cr
#'
#' Note: A wrong selection will causes a function error. Please change this
#' argument only if you have reasons to do so.
#' @param compatibility.mode \code{\link{logical}} (with default): this option
#' is parsed only if a BIN/BINX file is produced and it will reset all position
#' values to a max. value of 48, cf.\code{\link{write_R2BIN}}
#' @param txtProgressBar \code{\link{logical}} (with default): enables
#' \code{TRUE} or disables \code{FALSE} the progression bars during import and
#' export
#' @return An \code{\linkS4class{RLum.Results}} object is returned with the
#' following structure:\cr .. $irr.times (data.frame)\cr
#'
#' If a BINX-file path and name is set, the output will be additionally
#' transferred to a new BINX-file with the function name as suffix. For the
#' output the path of the input BINX-file itself is used. Note that this will
#' not work if the input object is a file path to an XSYG-file. In this case
#' the argument input is ignored.
#' @note The produced output object contains still the irradiation steps to
#' keep the output transparent. However, for the BINX-file export this steps
#' are removed as the BINX-file format description does not allow irradiations
#' as separat sequences steps.\cr
#'
#' Know issue: The 'fading correction' menu in the Analyst will not work appear
#' with the produced BIN/BINX-file due to hidden bits, which are not reproduced
#' by the function \code{write_R2BIN()} or if it appears it stops with a
#' floating point error. \cr
#'
#' Negative values for \code{TIMESINCELAS.STEP}? Yes, this is possible and no
#' bug, as in the XSYG file multiple curves are stored for one step. Example: A
#' TL step may comprise three curves: (a) counts vs. time, (b) measured
#' temperature vs. time and (c) predefined temperature vs. time. Three curves,
#' but they are all belonging to one TL measurement step, but with regard to
#' the time stamps this could produce negative values as the important function
#' (\code{\link{read_XSYG2R}}) do not change the order of entries for one step
#' towards a correct time order.
#' @section Function version: 0.2.1
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#' @seealso \code{\linkS4class{RLum.Analysis}},
#' \code{\linkS4class{RLum.Results}}, \code{\linkS4class{Risoe.BINfileData}},
#' \code{\link{read_XSYG2R}}, \code{\link{read_BIN2R}}, \code{\link{write_R2BIN}}
#' @references Duller, G., 2007. Analyst.
#' @keywords IO manip
#' @examples
#'
#'
#' ## (1) - example for your own data
#' ##
#' ## set files and run function
#' #
#' #   file.XSYG <- file.choose()
#' #   file.BINX <- file.choose()
#' #
#' #     output <- extract_IrradiationTimes(file.XSYG = file.XSYG, file.BINX = file.BINX)
#' #     get_RLum(output)
#' #
#' ## export results additionally to a CSV.file in the same directory as the XSYG-file
#' #       write.table(x = get_RLum(output),
#' #                   file = paste0(file.BINX,"_extract_IrradiationTimes.csv"),
#' #                   sep = ";",
#' #                   row.names = FALSE)
#'
#' @export
extract_IrradiationTimes <- function(
  object,
  file.BINX,
  recordType = c("irradiation (NA)", "IRSL (UVVIS)", "OSL (UVVIS)", "TL (UVVIS)"),
  compatibility.mode = TRUE,
  txtProgressBar = TRUE
){

  # Integrity tests -----------------------------------------------------------------------------

  ##check whether an character or an RLum.Analysis object is provided
  if(is(object)[1] != "character" & is(object)[1] != "RLum.Analysis"){

    stop("[extract_IrradiationTimes()] Input object is neither of type 'character' nor of type 'RLum.Analysis'.")

  }else if(is(object)[1] == "character"){

    ##set object to file.XSYG
    file.XSYG <- object

    ##XSYG
    ##check if file exists
    if(file.exists(file.XSYG) == FALSE){

      stop("[extract_IrradiationTimes()] Wrong XSYG file name or file does not exsits!")

    }

    ##check if file is XML file
    if(tail(unlist(strsplit(file.XSYG, split = "\\.")), 1) != "xsyg" &
         tail(unlist(strsplit(file.XSYG, split = "\\.")), 1) != "XSYG" ){

      stop("[extract_IrradiationTimes()] File is not of type 'XSYG'!")

    }

    ##BINX
    if(!missing(file.BINX)){

      ##check if file exists
      if(file.exists(file.BINX) == FALSE){

        stop("[extract_IrradiationTimes()] Wrong BINX file name or file does not exsits!")

      }

      ##check if file is XML file
      if(tail(unlist(strsplit(file.BINX, split = "\\.")), 1) != "binx" &
           tail(unlist(strsplit(file.BINX, split = "\\.")), 1) != "BINX" ){

        stop("[extract_IrradiationTimes()] File is not of type 'BINX'!")

      }

    }

    # Settings and import XSYG --------------------------------------------------------------------

    temp.XSYG <- read_XSYG2R(file.XSYG, txtProgressBar = txtProgressBar)

    if(!missing(file.BINX)){
      temp.BINX <- read_BIN2R(file.BINX, txtProgressBar = txtProgressBar)
      temp.BINX.dirname <- (dirname(file.XSYG))
    }


    # Some data preparation -----------------------------------------------------------------------
    ##set list
    temp.sequence.list <- list()

    ##select all analysis objects and combinde them
    for(i in 1:length(temp.XSYG)){

      ##select sequence and reduce the data set to really wanted values
      temp.sequence.list[[i]] <- get_RLum(temp.XSYG[[i]]$Sequence.Object,
                                                   recordType = recordType,
                                                   drop = FALSE)


      ##get corresponding position number, this will be needed later on
      temp.sequence.position <- as.numeric(as.character(temp.XSYG[[i]]$Sequence.Header["position",]))

    }



  }else{

    ##now we assume a single RLum.Analysis object
    ##select sequence and reduce the data set to really wanted values, note that no
    ##record selection was made!
    temp.sequence.list <- list(object)

  }





  ##merge objects
  if(length(temp.sequence.list)>1){

    temp.sequence <- merge_RLum(temp.sequence.list)

  }else{

    temp.sequence <- temp.sequence.list[[1]]

  }


  # Grep relevant information -------------------------------------------------------------------

  ##Sequence STEP
  STEP <- sapply(1:length_RLum(temp.sequence), function(x){

    get_RLum(temp.sequence, record.id = x)@recordType

  })

  #START time of each step
  temp.START <- unname(sapply(1:length_RLum(temp.sequence), function(x){

    get_RLum(get_RLum(temp.sequence, record.id = x), info.object = c("startDate"))

  }))

  ##DURATION of each STEP
  DURATION.STEP <- sapply(1:length_RLum(temp.sequence), function(x){

    max(get_RLum(get_RLum(temp.sequence, record.id = x))[,1])

  })



  ##a little bit reformatting.
  START <- strptime(temp.START, format = "%Y%m%d%H%M%S", tz = "GMT")

  ##Calculate END time of each STEP
  END <- START + DURATION.STEP

  ##add position number so far an XSYG file was the input
  if(exists("file.XSYG")){

    POSITION <- rep(temp.sequence.position, each = length_RLum(temp.sequence))

  }else if(!inherits(try(
    get_RLum(
      get_RLum(temp.sequence, record.id = 1), info.object = "position"),
    silent = TRUE), "try-error")){

    ##DURATION of each STEP
    POSITION <- unname(sapply(1:length_RLum(temp.sequence), function(x){

      get_RLum(get_RLum(temp.sequence, record.id = x),info.object = "position")

    }))

  }else{

    POSITION <- NA

  }


  ##Combine the results
  temp.results <- data.frame(POSITION,STEP,START,DURATION.STEP,END)


  # Calculate irradiation duration ------------------------------------------------------------

  ##set objects
  time.irr.duration <- NA

  IRR_TIME <- unlist(sapply(1:nrow(temp.results), function(x){

    if(temp.results[x,"STEP"] == "irradiation (NA)"){

      time.irr.duration <<- temp.results[x,"DURATION.STEP"]
      return(0)

    }else{

      if(is.na(time.irr.duration)){

        return(0)

      }else{

        return(time.irr.duration)

      }

    }

  }))


  # Calculate time since irradiation ------------------------------------------------------------

  ##set objects
  time.irr.end <- NA

  TIMESINCEIRR <- unlist(sapply(1:nrow(temp.results), function(x){

    if(temp.results[x,"STEP"] == "irradiation (NA)"){

      time.irr.end<<-temp.results[x,"END"]
      return(-1)

    }else{

      if(is.na(time.irr.end)){

        return(-1)

      }else{

        return(difftime(temp.results[x,"START"],time.irr.end, units = "secs"))

      }

    }

  }))



  # Calculate time since last step --------------------------------------------------------------


  TIMESINCELAST.STEP <- unlist(sapply(1:nrow(temp.results), function(x){

    if(x == 1){
      return(0)
    }else{
      return(difftime(temp.results[x,"START"],temp.results[x-1, "END"], units = "secs"))
    }


  }))


  # Combine final results -----------------------------------------------------------------------

  ##results table, export as CSV
  results <- cbind(temp.results,IRR_TIME, TIMESINCEIRR,TIMESINCELAST.STEP)

  # Write BINX-file if wanted -------------------------------------------------------------------
  if(!missing(file.BINX)){

    ##(1) remove all irradiation steps as there is no record in the BINX file and update information
    results.BINX <- results[-which(results[,"STEP"] == "irradiation (NA)"),]

    ##(1a)  update information
    temp.BINX@METADATA[,c("IRR_TIME", "TIMESINCEIRR")] <- results.BINX[,c("IRR_TIME","TIMESINCEIRR")]

    ##(2) compare entries in the BINX-file with the entries in the table to make sure
    ## that both have the same length
    if(!missing(file.BINX)){
      if(nrow(results.BINX) == nrow(temp.BINX@METADATA)){

        ##update BINX-file
        write_R2BIN(temp.BINX, version = "06",
                   file = paste0(file.BINX,"_extract_IrradiationTimes.BINX"),
                   compatibility.mode =  compatibility.mode,
                   txtProgressBar = txtProgressBar)


      }
    }else{

      warning("XSYG and BINX-file do not contain similar entries. BINX-file update skipped!")

    }
  }


  # Output --------------------------------------------------------------------------------------
  return(set_RLum(class = "RLum.Results", data = list(irr.times = results)))
}
