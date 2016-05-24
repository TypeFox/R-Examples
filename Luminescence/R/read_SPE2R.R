#' Import Princeton Intruments (TM) SPE-file into R
#'
#' Function imports Princeton Instruments (TM) SPE-files into R environment and
#' provides \code{RLum} objects as output.
#'
#' Function provides an import routine for the Princton Instruments SPE format.
#' Import functionality is based on the file format description provided by
#' Princton Instruments and a MatLab script written by Carl Hall (s.
#' references).
#'
#' @param file \link{character} (\bold{required}): spe-file name (including
#' path), e.g. \cr [WIN]: \code{read_SPE2R("C:/Desktop/test.spe")}, \cr
#' [MAC/LINUX]: \code{readSPER("/User/test/Desktop/test.spe")}
#'
#' @param output.object \code{\link{character}} (with default): set \code{RLum}
#' output object.  Allowed types are \code{"RLum.Data.Spectrum"},
#' \code{"RLum.Data.Image"} or \code{"matrix"}
#'
#' @param frame.range \code{\link{vector}} (optional): limit frame range, e.g.
#' select first 100 frames by \code{frame.range = c(1,100)}
#'
#' @param txtProgressBar \link{logical} (with default): enables or disables
#' \code{\link{txtProgressBar}}.
#'
#' @return Depending on the chosen option the functions returns three different
#' type of objects:\cr
#'
#' \code{output.object}. \cr
#'
#' \code{RLum.Data.Spectrum}\cr
#'
#' An object of type \code{\linkS4class{RLum.Data.Spectrum}} is returned.  Row
#' sums are used to integrate all counts over one channel.
#'
#' \code{RLum.Data.Image}\cr
#'
#' An object of type \code{\linkS4class{RLum.Data.Image}} is returned.  Due to
#' performace reasons the import is aborted for files containing more than 100
#' frames. This limitation can be overwritten manually by using the argument
#' \code{frame.frange}.
#'
#' \code{matrix}\cr
#'
#' Returns a matrix of the form: Rows = Channels, columns = Frames. For the
#' transformation the function \code{\link{get_RLum}} is used,
#' meaning that the same results can be obtained by using the function
#' \code{\link{get_RLum}} on an \code{RLum.Data.Spectrum} or \code{RLum.Data.Image} object.
#' @note \bold{The function does not test whether the input data are spectra or
#' pictures for spatial resolved analysis!}\cr
#'
#' The function has been successfully tested for SPE format versions 2.x.
#'
#' \emph{Currently not all information provided by the SPE format are
#' supported.}
#'
#' @section Function version: 0.1.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\link{readBin}}, \code{\linkS4class{RLum.Data.Spectrum}},
#' \code{\link[raster]{raster}}
#'
#' @references Princeton Instruments, 2014. Princeton Instruments SPE 3.0 File
#' Format Specification, Version 1.A,
#' \url{ftp://ftp.princetoninstruments.com/Public/Manuals/Princeton\%20Instruments/SPE\%203.0\%20File\%20Format\%20Specification.pdf}
#'
#' Hall, C., 2012: readSPE.m.
#' \url{http://www.mathworks.com/matlabcentral/fileexchange/35940-readspe/content/readSPE.m}
#'
#' @aliases readSPE2R
#'
#' @keywords IO
#'
#' @examples
#'
#'
#' ## to run examples uncomment lines and run the code
#'
#' ##(1) Import data as RLum.Data.Spectrum object
#' #file <- file.choose()
#' #temp <- read_SPE2R(file)
#' #temp
#'
#' ##(2) Import data as RLum.Data.Image object
#' #file <- file.choose()
#' #temp <- read_SPE2R(file, output.object = "RLum.Data.Image")
#' #temp
#'
#' ##(3) Import data as matrix object
#' #file <- file.choose()
#' #temp <- read_SPE2R(file, output.object = "matrix")
#' #temp
#'
#' ##(4) Export raw data to csv, if temp is a RLum.Data.Spectrum object
#' # write.table(x = get_RLum(temp),
#' #             file = "[your path and filename]",
#' #             sep = ";", row.names = FALSE)
#'
#'
#' @export
read_SPE2R <- function(
  file,
  output.object = "RLum.Data.Image",
  frame.range,
  txtProgressBar = TRUE
){

  # Consistency check -------------------------------------------------------

  ##check if file exists
  if(file.exists(file) == FALSE){

    stop("[read_SPE2R()] File not found!")

  }

  ##check file extension
  if(strsplit(file, split = "\\.")[[1]][2] != "SPE"){

    temp.text <- paste("[read_SPE2R()] Unsupported file format: *.",
                       strsplit(file, split = "\\.")[[1]][2], sep = "")

    stop(temp.text)

  }


  # Open Connection ---------------------------------------------------------

  #open connection
  con<-file(file, "rb")

  # read header -------------------------------------------------------------

  temp <- readBin(con, what="int", 2, size=2, endian="little", signed = TRUE)
  ControllerVersion <- temp[1] #Hardware version
  LogicOutput <- temp[2] #Definition of Output BNC

  temp <- readBin(con, what="int", 2, size=2, endian="little", signed = FALSE)
  AmpHiCapLowNoise <- temp[1] #Amp Switching Mode
  xDimDet <- temp[2] #Detector x dimension of chip.

  #timing mode
  mode <- readBin(con, what="int", 1, size=2, endian="little", signed = TRUE)

  #alternative exposure, in sec.
  exp_sec <- readBin(con, what="double", 1, size=4, endian="little")

  temp <- readBin(con, what="int", 2, size=2, endian="little", signed = TRUE)
  VChipXdim <- temp[1] # Virtual Chip X dim
  VChipYdim <- temp[2] # Virtual Chip Y dim

  #y dimension of CCD or detector.
  yDimDet <- readBin(con, what="int", 1, size=2, endian="little", signed = TRUE)

  #Date
  Date <- readChar(con, 10, useBytes=TRUE)

  ##jump
  stepping <- readBin(con, what="raw", 4, size=1, endian="little", signed = TRUE)

  #Old number of scans - should always be -1
  noscan <- readBin(con, what="int", 1, size=2, endian="little", signed = TRUE)

  #Detector Temperature Set
  DetTemperature <- readBin(con, what="double", 1, size=4, endian="little")

  # CCD/DiodeArray type
  DetType <- readBin(con, what="int", 1, size=2, endian="little", signed = TRUE)

  #actual # of pixels on x axis
  xdim <- readBin(con, what="int", 1, size=2, endian="little", signed = FALSE)

  ##jump
  stepping <- readBin(con, what="raw", 64, size=1, endian="little", signed = TRUE)

  ##experiment data type
  ##0 = 32f (4 bytes)
  ##1 = 32s (4 bytes)
  ##3 = 16u (2 bytes)
  ##8 = 32u (4 bytes)
  datatype <- readBin(con, what="int", 1, size=2, endian="little", signed = TRUE)

  ##jump
  stepping <- readBin(con, what="raw", 546, size=1, endian="little")

  #y dimension of raw data.
  ydim <- readBin(con, what="int", 1, size=2, endian="little", signed = FALSE)

  ##0=scrambled,1=unscrambled
  scramble <- readBin(con, what="int", 1, size=2, endian="little", signed = FALSE)

  ##jump
  stepping <- readBin(con, what="raw", 4, size=1, endian="little")

  #Number of scans (Early WinX)
  lnoscan <- readBin(con, what="int", 1, size=4, endian="little", signed = TRUE)

  #Number of Accumulations
  lavgexp <- readBin(con, what="int", 1, size=4, endian="little", signed = TRUE)

  ##Experiment readout time
  ReadoutTime <- readBin(con, what="double", 1, size=4, endian="little")

  #T/F Triggered Timing Option
  TriggeredModeFlag <- readBin(con, what="int", 1, size=2, endian="little", signed = TRUE)

  ##jump
  stepping <- readBin(con, what="raw", 768, size=1, endian="little")

  ##number of frames in file.
  NumFrames <- readBin(con, what="int", 1, size=4, endian="little", signed = TRUE)

  if(NumFrames > 100 & missing(frame.range) & output.object == "RLum.Data.Image"){

    error.message <- paste0("[read_SPE2R()] Import aborted. This file containes > 100 (", NumFrames, "). Use argument 'frame.range' to force import.")
    stop(error.message)

  }

  ##set frame.range
  if(missing(frame.range) == TRUE){frame.range <- c(1,NumFrames)}

  ##jump
  stepping <- readBin(con, what="raw", 542, size=1, endian="little")

  #file_header_ver
  file_header_ver <- readBin(con, what="double", 1, size=4, endian="little")

  ##jump
  stepping <- readBin(con, what="raw", 1000, size=1, endian="little")

  ##WinView_id - set to 19,088,743 (or 1234567 hex) (required for legacy reasons)
  WinView_id <- readBin(con, what="integer", 1, size=4, endian="little", signed = TRUE)

  ##jump
  stepping <- readBin(con, what="raw", 1098, size=1, endian="little")

  ##lastvalue - set to 21,845 (or 5555 hex) (required for legacy reasons)
  lastvalue <- readBin(con, what="integer", 1, size=2, endian="little", signed = TRUE)


  ##end header
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ##create info element list from data
  temp.info <- list(ControllerVersion,
                    LogicOutput,
                    AmpHiCapLowNoise,
                    xDimDet, yDimDet,
                    xdim, ydim,
                    VChipXdim, VChipYdim,
                    Date,
                    noscan,
                    mode,  exp_sec,
                    DetTemperature,
                    DetType,
                    datatype,
                    scramble,
                    lnoscan,
                    lavgexp,
                    ReadoutTime,
                    TriggeredModeFlag,
                    NumFrames,
                    file_header_ver)

  ##set name for list elements
  names(temp.info) <- c("ControllerVersion", "LogicOutput", "AmpHiCapLowNoise", "xDimDet", "yDimDet",
                        "xdim", "ydim", "VChipXdim", "VChipYdim", "Date", "noscan", "mode", "exp_sec",
                        "DetTemperature", "DetType", "datatype", "scramble", "lnoscan", "lavgexp",
                        "ReadoutTime", "TriggeredModeFlag", "NumFrames", "file_header_ver")

  # read count value data ---------------------------------------------------
  ##set functions

  if(datatype  == 0){

    read.data <- function(n.counts){
      readBin(con, what="double", n.counts, size=4, endian="little")
    }

  }else if(datatype == 1){

    read.data <- function(n.counts){
      readBin(con, what="integer", n.counts, size=4, endian="little", signed = TRUE)
    }

  }else if(datatype == 2){

    read.data <- function(n.counts){
      readBin(con, what="integer", n.counts, size=2, endian="little", signed = TRUE)
    }

  }else if(datatype == 3){

    read.data <- function(n.counts){

      readBin(con, what="int", n.counts, size=2, endian="little", signed = FALSE)

    }

  }else if(datatype == 8){

    read.data <- function(n.counts){
      readBin(con, what="integer", n.counts, size=4, endian="little", signed = FALSE)
    }

  }else{

    stop("[read_SPE2R()] Unknown 'datatype'.")

  }


  ##loop over all frames
  ##output
  cat(paste("\n[read_SPE2R.R]\n\t >> ",file,sep=""), fill=TRUE)

  ##set progressbar
  if(txtProgressBar==TRUE){
    pb<-txtProgressBar(min=0,max=diff(frame.range)+1, char="=", style=3)
  }

  ##stepping for frame range
  temp <- readBin(con, what = "raw", (min(frame.range)-1)*2, size = 1, endian = "little")

  for(i in 1:(diff(frame.range)+1)){#NumFrames

    temp.data <- matrix(read.data(n.counts = (xdim * ydim)),
                        ncol = ydim,
                        nrow = xdim)

    if(exists("data.list") == FALSE){

      data.list <- list(temp.data)

    }else{

      data.list <- c(data.list, list(temp.data))

    }

    ##update progress bar
    if(txtProgressBar==TRUE){
      setTxtProgressBar(pb, i)
    }

  }

  ##close
  if(txtProgressBar==TRUE){close(pb)

                           ##output
                           cat(paste("\t >> ",i," records have been read successfully!\n\n", sep=""))
  }

  # Output ------------------------------------------------------------------

  if(output.object == "RLum.Data.Spectrum" | output.object == "matrix"){

    ##to create a spectrum object the matrix has to transposed and
    ##the row sums are needed

    data.spectrum.vector <- sapply(1:length(data.list), function(x){

      rowSums(data.list[[x]])

    })

    ##split vector to matrix
    data.spectrum.matrix <- matrix(data.spectrum.vector,
                                   nrow = xdim,
                                   ncol = length(data.list))

    ##set column and row names
    colnames(data.spectrum.matrix) <- as.character(1:ncol(data.spectrum.matrix))
    rownames(data.spectrum.matrix) <- as.character(1:nrow(data.spectrum.matrix))


    ##set output object
    object <- set_RLum(
      class = "RLum.Data.Spectrum",
      originator = "read_SPE2R",
      recordType = "Spectrum",
                                     curveType = "measured",
                                     data = data.spectrum.matrix,
                                     info = temp.info)

    ##optional matrix object
    if(output.object == "matrix"){

      object <- get_RLum(object)}


  }else if(output.object == "RLum.Data.Image"){

    ##combine to raster
    data.raster.list <- lapply(1:length(data.list), function(x){

      if(txtProgressBar==TRUE){

        cat(paste("\r Converting to RasterLayer: ", x, "/",length(data.list), sep = ""))

      }

      raster::raster(t(data.list[[x]]),
             xmn = 0, xmx = max(xdim),
             ymn = 0, ymx = max(ydim))


    })

    ##Convert to raster brick
    data.raster <- raster::brick(x = data.raster.list)

    ##Create RLum.object
    object <- set_RLum(
      class = "RLum.Data.Image",
      originator = "read_SPE2R",
      recordType = "Image",
                                  curveType = "measured",
                                  data = data.raster,
                                  info = temp.info)

  }else{

    stop("[read_SPE2R()] Chosen 'output.object' not supported. Please check manual!")

  }

  ##close con
  close(con)


  ##return values
  return(object)

}

## ---- DEPRECATED GENERICS
# .Deprecated in package version 0.5.0
# .Defunct in 0.5.1
# Removed in 0.6.0
#' @noRd
#' @export
readSPE2R <- function(...) {
  .Defunct("read_SPE2R")
  read_SPE2R(...)
}


