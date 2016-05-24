#' Import Risoe BIN-file into R
#'
#' Import a *.bin or a *.binx file produced by a Risoe DA15 and DA20 TL/OSL
#' reader into R.
#'
#' The binary data file is parsed byte by byte following the data structure
#' published in the Appendices of the Analyst manual p. 42.\cr\cr For the
#' general BIN-file structure, the reader is referred to the Risoe website:
#' \code{http://www.nutech.dtu.dk/}
#'
#' @param file \code{\link{character}} or \code{\link{list}} (\bold{required}): path and file name of the
#' BIN/BINX file. If input is a \code{list} it should comprise only \code{character}s representing
#' each valid path and BIN/BINX-file names.
#' Alternatively the input character can be just a directory (path), in this case the
#' the function tries to detect and import all BIN/BINX files found in the directory.
#'
#' @param show.raw.values \link{logical} (with default): shows raw values from
#' BIN file for \code{LTYPE}, \code{DTYPE} and \code{LIGHTSOURCE} without
#' translation in characters.
#'
#' @param n.records \link{raw} (optional): limits the number of imported
#' records. Can be used in combination with \code{show.record.number} for
#' debugging purposes, e.g. corrupt BIN files.
#'
#' @param duplicated.rm \code{\link{logical}} (with default): remove duplicated entries if \code{TRUE}.
#' This may happen due to an erroneous produced BIN/BINX-file. This option compares only
#' predeccessor and successor.
#'
#' @param position \code{\link{numeric}} (optional): imports only the selected position. Note:
#' the import performance will not benefit by any selection made here.
#'
#' @param fastForward \code{\link{logical}} (with default): if \code{TRUE} for a
#' more efficient data processing only a list of \code{RLum.Analysis} objects is returned instead
#' of a \link{Risoe.BINfileData-class} object
#'
#' @param show.record.number \link{logical} (with default): shows record number
#' of the imported record, for debugging usage only.
#'
#' @param txtProgressBar \link{logical} (with default): enables or disables
#' \code{\link{txtProgressBar}}.
#'
#' @param forced.VersionNumber \link{integer} (optional): allows to cheat the
#' version number check in the function by own values for cases where the
#' BIN-file version is not supported.\cr Note: The usage is at own risk, only
#' supported BIN-file versions have been tested.
#'
#' @param \dots further arguments that will be passed to the function
#' \code{\link{Risoe.BINfileData2RLum.Analysis}}. Please note that any matching argument
#' automatically sets \code{fastForward = TRUE}
#'
#' @return Returns an S4 \link{Risoe.BINfileData-class} object containing two
#' slots:\cr \item{METADATA}{A \link{data.frame} containing all variables
#' stored in the bin-file.} \item{DATA}{A \link{list} containing a numeric
#' \link{vector} of the measured data. The ID corresponds to the record ID in
#' METADATA.}\cr
#'
#' If \code{fastForward = TRUE} a list of \code{\linkS4class{RLum.Analysis}} object is returned. The
#' internal coercing is done using the function \code{\link{Risoe.BINfileData2RLum.Analysis}}
#'
#'
#' @note The function works for BIN/BINX-format versions 03, 04, 06 and 07. The
#' version number depends on the used Sequence Editor.\cr\cr \bold{Potential
#' other BIN/BINX-format versions are currently not supported. The
#' implementation of version 07 support could not been tested so far.}.
#'
#'
#' @section Function version: 0.12.0
#'
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France), Margret C. Fuchs, HZDR Freiberg, (Germany)
#'
#'
#' @seealso \code{\link{write_R2BIN}}, \code{\linkS4class{Risoe.BINfileData}},
#' \code{\link[base]{readBin}}, \code{\link{merge_Risoe.BINfileData}}, \code{\linkS4class{RLum.Analysis}}
#' \code{\link[utils]{txtProgressBar}}
#'
#'
#' @references Duller, G., 2007. Analyst.
#' \url{http://www.nutech.dtu.dk/english/~/media/Andre_Universitetsenheder/Nutech/Produkter\%20og\%20services/Dosimetri/radiation_measurement_instruments/tl_osl_reader/Manuals/analyst_manual_v3_22b.ashx}
#'
#' @aliases readBIN2R
#'
#' @keywords IO
#'
#'
#' @examples
#'
#'
#' ##(1) import Risoe BIN-file to R (uncomment for usage)
#'
#' #FILE <- file.choose()
#' #temp <- read_BIN2R(FILE)
#' #temp
#'
#' @export
read_BIN2R <- function(
  file,
  show.raw.values = FALSE,
  position = NULL,
  n.records = NULL,
  duplicated.rm = FALSE,
  fastForward = FALSE,
  show.record.number = FALSE,
  txtProgressBar = TRUE,
  forced.VersionNumber = NULL,
  ...
){

  # Self Call -----------------------------------------------------------------------------------
  # Option (a): Input is a list, every element in the list will be treated as file connection
  # with that many file can be read in at the same time
  # Option (b): The input is just a path, the function tries to grep ALL BIN/BINX files in the
  # directory and import them, if this is detected, we proceed as list

  if(is(file, "character")) {

    ##If this is not really a path we skip this here
    if (dir.exists(file) & length(dir(file)) > 0) {
      cat("[read_BIN2R()] Directory detected, trying to extract '*.bin'/'*.binx' files ...\n")
      file <-
        as.list(c(
          paste0(file,dir(
            file, recursive = FALSE, pattern = ".bin"
          )),
          paste0(file,dir(
            file, recursive = FALSE, pattern = ".binx"
          )),
          paste0(file,dir(
            file, recursive = FALSE, pattern = ".BIN"
          )),  paste0(file,dir(
            file, recursive = FALSE, pattern = ".BINX"
          ))
        ))

    }

  }

  if (is(file, "list")) {
    temp.return <- lapply(1:length(file), function(x) {
      temp <- read_BIN2R(
        file = file[[x]],
        fastForward = fastForward,
        position = position,
        n.records = n.records,
        show.raw.values =  show.raw.values,
        show.record.number = show.record.number,
        txtProgressBar = txtProgressBar,
        forced.VersionNumber = forced.VersionNumber,
        ...
      )

    })

    ##return
    if (fastForward) {
      return(unlist(temp.return, recursive = FALSE))

    }else{
      return(temp.return)

    }

  }



  # Integrity checks ------------------------------------------------------

  ##check if file exists
  if(!file.exists(file)){

    stop("[read_BIN2R()] File does not exists!")

  }

  ##check if file is a BIN or BINX file
  if(!(TRUE%in%(c("BIN", "BINX", "bin", "binx")%in%tail(
    unlist(strsplit(file, split = "\\.")), n = 1)))){

    stop("[read_BIN2R()] Input is not a file or not of type 'BIN' or 'BINX'!")

  }

  # Config ------------------------------------------------------------------

  ##set supported BIN format version
  VERSION.supported <- as.raw(c(03, 04, 06, 07))


  # Short file parsing to get number of records -------------------------------------------------

  #open connection
  con<-file(file, "rb")

  ##get information about file size
  file.size<-file.info(file)

  ##read data up to the end of con

  ##set ID
  temp.ID<-0

  ##start for BIN-file check up
  while(length(temp.VERSION<-readBin(con, what="raw", 1, size=1, endian="litte"))>0) {

    ##force version number
    if(!is.null(forced.VersionNumber)){
      temp.VERSION <- as.raw(forced.VersionNumber)
    }

    ##stop input if wrong VERSION
    if((temp.VERSION%in%VERSION.supported) == FALSE){

      ##close connection
      close(con)

      ##show error message
      error.text <- paste("[read_BIN2R()] The BIN-format version (",temp.VERSION,") of this file is currently not supported! Supported version numbers are: ",paste(VERSION.supported,collapse=", "),".",sep="")

      stop(error.text)

    }

    #empty byte position
    EMPTY<-readBin(con, what="raw", 1, size=1, endian="litte")

    if(temp.VERSION==06 | temp.VERSION==07){

      ##GET record LENGTH
      temp.LENGTH  <- readBin(con, what="int", 1, size=4, endian="little")
      STEPPING <- readBin(con, what="raw", temp.LENGTH-6, size=1, endian="litte")

    }else{

      ##GET record LENGTH
      temp.LENGTH  <- readBin(con, what="int", 1, size=2, endian="little")
      STEPPING <- readBin(con, what="raw", temp.LENGTH-4, size=1, endian="litte")

    }

    temp.ID<-temp.ID+1

  }

  ##close con
  close(con)

# Set Lookup tables  --------------------------------------------------------------------------

  ##LTYPE
  LTYPE.lookup <- c(
    "0" = "TL",
    "1" = "OSL",
    "2" = "IRSL",
    "3" = "M-IR",
    "4" = "M-VIS",
    "5" = "TOL",
    "6" = "TRPOSL",
    "7" = "RIR",
    "8" = "RBR",
    "9" = "USER",
    "10" = "POSL",
    "11" = "SGOSL",
    "12" = "RL",
    "13" = "XRF"
  )

  ##DTYPE
  DTYPE.lookup <-
    c(
      "0" = "Natural",
      "1" = "N+dose",
      "2" = "Bleach",
      "3" = "Bleach+dose",
      "4" = "Natural (Bleach)",
      "5" = "N+dose (Bleach)",
      "6" = "Dose",
      "7" = "Background"
    )

  ##LIGHTSOURCE
  LIGHTSOURCE.lookup <- c(
    "0" = "None",
    "1" = "Lamp",
    "2" = "IR diodes/IR Laser",
    "3" = "Calibration LED",
    "4" = "Blue Diodes",
    "5" = "White light",
    "6" = "Green laser (single grain)",
    "7" = "IR laser (single grain)"
  )


  ##PRESET VALUES
  temp.CURVENO <- NA
  temp.FNAME <- NA
  temp.MEASTEMP <- NA
  temp.IRR_UNIT <- NA
  temp.IRR_DOSERATE <- NA
  temp.IRR_DOSERATEERR <- NA
  temp.TIMESINCEIRR <- NA
  temp.TIMETICK <- NA
  temp.ONTIME <- NA
  temp.OFFTIME <- NA
  temp.STIMPERIOD <- NA
  temp.GATE_ENABLED <- raw(length = 1)
  temp.ENABLE_FLAGS <- raw(length = 1)
  temp.GATE_START <- NA
  temp.GATE_STOP <- NA
  temp.GATE_END <- NA
  temp.PTENABLED <- raw(length = 1)
  temp.DTENABLED <- raw(length = 1)
  temp.DEADTIME <- NA
  temp.MAXLPOWER <- NA
  temp.XRF_ACQTIME <- NA
  temp.XRF_HV <- NA
  temp.XRF_CURR <- NA
  temp.XRF_DEADTIMEF <- NA
  temp.DETECTOR_ID <- NA
  temp.LOWERFILTER_ID <- NA
  temp.UPPERFILTER_ID <- NA
  temp.ENOISEFACTOR <- NA
  temp.SEQUENCE <- NA

  ##SET length of entire record
  n.length = temp.ID
  rm(temp.ID)

  ##initialise data.frame
  results.METADATA <- data.table(

    ID = integer(length = n.length),
    SEL = logical(length = n.length),
    VERSION = numeric(length = n.length),
    LENGTH = integer(length = n.length),
    PREVIOUS = integer(length = n.length),
    NPOINTS = integer(length = n.length),

    RUN = integer(length = n.length),
    SET = integer(length = n.length),
    POSITION = integer(length = n.length),
    GRAIN = integer(length = n.length),
    GRAINNUMBER = integer(length = n.length),
    CURVENO = integer(length = n.length),
    XCOORD = integer(length = n.length),
    YCOORD = integer(length = n.length),
    SAMPLE = character(length = n.length),
    COMMENT = character(length = n.length),

    SYSTEMID = integer(length = n.length),
    FNAME = character(length = n.length),
    USER = character(length = n.length),
    TIME = character(length = n.length),
    DATE = character(length = n.length),

    DTYPE = character(length = n.length),
    BL_TIME = numeric(length = n.length),
    BL_UNIT = integer(length = n.length),
    NORM1 = numeric(length = n.length),
    NORM2 = numeric(length = n.length),
    NORM3 = numeric(length = n.length),
    BG = numeric(length = n.length),
    SHIFT = integer(length = n.length),
    TAG = integer(length = n.length),

    LTYPE = character(length = n.length),
    LIGHTSOURCE = character(length = n.length),
    LPOWER = numeric(length = n.length),
    LIGHTPOWER = numeric(length = n.length),
    LOW = numeric(length = n.length),
    HIGH = numeric(length = n.length),
    RATE = numeric(length = n.length),
    TEMPERATURE = numeric(length = n.length),
    MEASTEMP = numeric(length = n.length),
    AN_TEMP = numeric(length = n.length),
    AN_TIME = numeric(length = n.length),
    TOLDELAY = integer(length = n.length),
    TOLON = integer(length = n.length),
    TOLOFF = integer(length = n.length),
    IRR_TIME = numeric(length = n.length),
    IRR_TYPE = integer(length = n.length),
    IRR_UNIT = integer(length = n.length),
    IRR_DOSERATE = numeric(length = n.length),
    IRR_DOSERATEERR = numeric(length = n.length),
    TIMESINCEIRR = numeric(length = n.length),
    TIMETICK = numeric(length = n.length),
    ONTIME = numeric(length = n.length),
    OFFTIME = numeric(length = n.length),
    STIMPERIOD = integer(length = n.length),
    GATE_ENABLED = numeric(length = n.length),
    ENABLE_FLAGS = numeric(length = n.length),
    GATE_START  = numeric(length = n.length),
    GATE_STOP = numeric(length = n.length),
    PTENABLED = numeric(length = n.length),
    DTENABLED = numeric(length = n.length),
    DEADTIME = numeric(length = n.length),
    MAXLPOWER = numeric(length = n.length),
    XRF_ACQTIME = numeric(length = n.length),
    XRF_HV = numeric(length = n.length),
    XRF_CURR = numeric(length = n.length),
    XRF_DEADTIMEF = numeric(length = n.length),

    DETECTOR_ID = integer(length = n.length),
    LOWERFILTER_ID = integer(length = n.length),
    UPPERFILTER_ID = integer(length = n.length),
    ENOISEFACTOR = numeric(length = n.length),

    SEQUENCE = character(length = n.length)

  ) #end set data table

  #set variable for DPOINTS handling
  results.DATA<-list()

  ##set list for RESERVED values
  results.RESERVED <- rep(list(list()), n.length)

  # Open Connection ---------------------------------------------------------

  ##show warning if version number check has been cheated

  if(!is.null(forced.VersionNumber)){
    warning("Argument 'forced.VersionNumber' has been used. BIN-file version might be not supported!")
  }

  #open connection
  con<-file(file, "rb")

  ##get information about file size
  file.size<-file.info(file)

  ##output
  cat(paste("\n[read_BIN2R()]\n\t >> ",file,sep=""), fill=TRUE)

  ##set progressbar
  if(txtProgressBar==TRUE){
    pb<-txtProgressBar(min=0,max=file.size$size, char="=", style=3)
  }

  ##read data up to the end of con

  ##set ID
  temp.ID<-0


  # LOOP --------------------------------------------------------------------

  ##start loop for import BIN data
  while(length(temp.VERSION<-readBin(con, what="raw", 1, size=1, endian="litte"))>0) {

    ##force version number
    if(!is.null(forced.VersionNumber)){
      temp.VERSION <- as.raw(forced.VersionNumber)
    }

    ##stop input if wrong VERSION
    if((temp.VERSION%in%VERSION.supported) == FALSE){

      ##close connection
      close(con)

      ##show error message
      error.text <- paste("[read_BIN2R()] The BIN-format version (",temp.VERSION,") of this file is currently not supported! Supported version numbers are: ",paste(VERSION.supported,collapse=", "),".",sep="")

      stop(error.text)

    }

    ##print record ID for debugging purposes
    if(show.record.number == TRUE){



      cat(temp.ID,",", sep = "")
      if(temp.ID%%10==0){
        cat("\n")
      }
    }


    #empty byte position
    EMPTY<-readBin(con, what="raw", 1, size=1, endian="litte")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # BINX FORMAT SUPPORT -----------------------------------------------------
    if(temp.VERSION==06 | temp.VERSION==07){

      ##(1) Header size and strucutre
      ##LENGTH, PREVIOUS, NPOINTS, LTYPE
      temp <- readBin(con, what="int", 3, size=4, endian="little")

      temp.LENGTH <- temp[1]
      temp.PREVIOUS <- temp[2]
      temp.NPOINTS <- temp[3]

      ##(2) Sample characteristics
      ##RUN, SET, POSITION, GRAINNUMBER, CURVENO, XCOORD, YCOORD
      temp <- readBin(con, what="int", 7, size=2, endian="little")

      temp.RUN <- temp[1]
      temp.SET <- temp[2]
      temp.POSITION <- temp[3]
      temp.GRAINNUMBER <- temp[4]
      temp.CURVENO <- temp[5]
      temp.XCOORD <- temp[6]
      temp.YCOORD <- temp[7]

      ##SAMPLE, COMMENT
      ##SAMPLE
      SAMPLE_SIZE<-readBin(con, what="int", 1, size=1, endian="little")
      temp.SAMPLE<-readChar(con, SAMPLE_SIZE, useBytes=TRUE)
      #however it should be set to 20

      #step forward in con
      if(20-c(SAMPLE_SIZE)>0){
        STEPPING<-readBin(con, what="raw", (20-c(SAMPLE_SIZE)),
                          size=1, endian="little")
      }

      ##COMMENT
      COMMENT_SIZE<-readBin(con, what="int", 1, size=1, endian="little")
      temp.COMMENT<-readChar(con, COMMENT_SIZE, useBytes=TRUE) #set to 80 (manual)

      #step forward in con
      if(80-c(COMMENT_SIZE)>0){
        STEPPING<-readBin(con, what="raw", (80-c(COMMENT_SIZE)),
                          size=1, endian="little")
      }

      ##(3) Instrument and sequence characteristic
      ##SYSTEMID
      temp.SYSTEMID <- readBin(con, what="int", 1, size=2, endian="little")

      ##FNAME
      FNAME_SIZE<-readBin(con, what="int", 1, size=1, endian="little")

      ##correct for 0 file name length
      if(length(FNAME_SIZE)>0){
        temp.FNAME<-readChar(con, FNAME_SIZE, useBytes=TRUE) #set to 100 (manual)
      }else{
        FNAME_SIZE <- 0
      }

      #step forward in con
      if(100-c(FNAME_SIZE)>0){
        STEPPING<-readBin(con, what="raw", (100-c(FNAME_SIZE)),
                          size=1, endian="little")
      }

      ##USER
      USER_SIZE<-readBin(con, what="int", 1, size=1, endian="little")

      ##correct for 0 user size length
      if (length(USER_SIZE) > 0) {
        temp.USER <-
          readChar(con, USER_SIZE, useBytes = TRUE) #set to 30 (manual)
      }else{
        USER_SIZE <- 0

      }

      #step forward in con
      if(30-c(USER_SIZE)>0){
        STEPPING<-readBin(con, what="raw", (30-c(USER_SIZE)),
                          size=1, endian="little")
      }

      ##TIME
      TIME_SIZE<-readBin(con, what="int", 1, size=1, endian="little")

      ##time size corrections for wrong time formats; set n to 6 for all values
      ##accoording the handbook of Geoff Duller, 2007
      if(length(TIME_SIZE)>0){
        temp.TIME<-readChar(con, TIME_SIZE, useBytes=TRUE)
      }else{
        TIME_SIZE <- 0

      }

      if(6-TIME_SIZE>0){

        STEPPING<-readBin(con, what="raw", (6-TIME_SIZE),
                          size=1, endian="little")
      }


      ##DATE
      DATE_SIZE<-readBin(con, what="int", 1, size=1, endian="little")

      ##date size corrections for wrong date formats; set n to 6 for all values
      ##accoording the handbook of Geoff Duller, 2007
      DATE_SIZE<-6
      temp.DATE<-readChar(con, DATE_SIZE, useBytes=TRUE)


      ##(4) Analysis

      ##DTYPE
      temp.DTYPE<-readBin(con, what="int", 1, size=1, endian="little")

      ##BL_TIME
      temp.BL_TIME<-readBin(con, what="double", 1, size=4, endian="little")

      ##BL_UNIT
      temp.BL_UNIT<-readBin(con, what="int", 1, size=1, endian="little")

      ##NORM1, NORM2, NORM3, BG
      temp <- readBin(con, what="double", 4, size=4, endian="little")

      temp.NORM1 <- temp[1]
      temp.NORM2 <- temp[2]
      temp.NORM3 <- temp[3]
      temp.BG <- temp[4]

      ##SHIFT
      temp.SHIFT<- readBin(con, what="integer", 1, size=2, endian="little")

      ##TAG
      temp.TAG <- readBin(con, what="int", 1, size=1, endian="little")

      ##RESERVED
      temp.RESERVED1 <-readBin(con, what="raw", 20, size=1, endian="little")

      ##(5) Measurement characteristics

      ##LTYPE
      temp.LTYPE <- readBin(con, what="int", 1, size=1, endian="little")

      ##LTYPESOURCE
      temp.LIGHTSOURCE <- readBin(con, what="int", 1, size=1, endian="little")

      ##LIGHTPOWER, LOW, HIGH, RATE
      temp <- readBin(con, what="double", 4, size=4, endian="little")

      temp.LIGHTPOWER <- temp[1]
      temp.LOW <- temp[2]
      temp.HIGH <- temp[3]
      temp.RATE <- temp[4]

      ##TEMPERATURE
      temp.TEMPERATURE <- readBin(con, what="int", 1, size=2, endian="little")

      ##MEASTEMP
      temp.MEASTEMP <- readBin(con, what="integer", 1, size=2, endian="little")

      ##AN_TEMP
      temp.AN_TEMP <- readBin(con, what="double", 1, size=4, endian="little")

      ##AN_TIME
      temp.AN_TIME <- readBin(con, what="double", 1, size=4, endian="little")

      ##DELAY, ON, OFF
      temp <- readBin(con, what="int", 3, size=2, endian="little")

      temp.TOLDELAY <- temp[1]
      temp.TOLON <- temp[2]
      temp.TOLOFF <- temp[3]

      ##IRR_TIME
      temp.IRR_TIME <- readBin(con, what="double", 1, size=4, endian="little")

      ##IRR_TYPE
      temp.IRR_TYPE <- readBin(con, what="int", 1, size=1, endian="little")

      ##IRR_DOSERATE
      temp.IRR_DOSERATE <- readBin(con, what="double", 1, size=4, endian="little")

      ##IRR_DOSERATEERR
      temp.IRR_DOSERATEERR <- readBin(con, what="double", 1, size=4, endian="little")

      ##TIMESINCEIRR
      temp.TIMESINCEIRR <- readBin(con, what="integer", 1, size=4, endian="little")

      ##TIMETICK
      temp.TIMETICK <- readBin(con, what="double", 1, size=4, endian="little")

      ##ONTIME
      temp.ONTIME <- readBin(con, what="integer", 1, size=4, endian="little")

      ##STIMPERIOD
      temp.STIMPERIOD <- readBin(con, what="integer", 1, size=4, endian="little")

      ##GATE_ENABLED
      temp.GATE_ENABLED <- readBin(con, what="raw", 1, size=1, endian="little")

      ##GATE_START
      temp.GATE_START <- readBin(con, what="integer", 1, size=4, endian="little")

      ##GATE_STOP
      temp.GATE_STOP <- readBin(con, what="integer", 1, size=4, endian="little")

      ##PTENABLED
      temp.PTENABLED <- readBin(con, what="raw", 1, size=1, endian="little")

      ##DTENABLED
      temp.DTENABLED <- readBin(con, what="raw", 1, size=1, endian="little")

      ##DEADTIME, MAXLPOWER, XRF_ACQTIME, XRF_HV
      temp <- readBin(con, what="double", 4, size=4, endian="little")

      temp.DEADTIME <- temp[1]
      temp.MAXLPOWER <- temp[2]
      temp.XRF_ACQTIME <- temp[3]
      temp.XRF_HV <- temp[4]

      ##XRF_CURR
      temp.XRF_CURR <- readBin(con, what="integer", 1, size=4, endian="little")

      ##XRF_DEADTIMEF
      temp.XRF_DEADTIMEF <- readBin(con, what="double", 1, size=4, endian="little")

      ###Account for differences between V6 and V7
      if(temp.VERSION == 06){

        ##RESERVED
        temp.RESERVED2<-readBin(con, what="raw", 24, size=1, endian="little")


      }else{

        ##DETECTOR_ID
        temp.DETECTOR_ID <- readBin(con, what="int", 1, size=1, endian="little")

        ##LOWERFILTER_ID, UPPERFILTER_ID
        temp <- readBin(con, what="int", 2, size=2, endian="little")

        temp.LOWERFILTER_ID <- temp[1]
        temp.UPPERFILTER_ID <- temp[2]

        ##ENOISEFACTOR
        temp.ENOISEFACTOR <- readBin(con, what="double", 1, size=4, endian="little")

        ##RESERVED
        temp.RESERVED2<-readBin(con, what="raw", 15, size=1, endian="little")

      }

      #DPOINTS
      temp.DPOINTS<-readBin(con, what="integer", temp.NPOINTS, size=4, endian="little")

    }else if(temp.VERSION == 04 | temp.VERSION == 03){
      ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ##START BIN FILE FORMAT SUPPORT  (vers. 03 and 04)
      ##LENGTH, PREVIOUS, NPOINTS, LTYPE

      temp <- readBin(con, what="int", 3, size=2, endian="little")

      temp.LENGTH <- temp[1]
      temp.PREVIOUS <- temp[2]
      temp.NPOINTS <- temp[3]

      ##LTYPE
      temp.LTYPE<-readBin(con, what="int", 1, size=1, endian="little")


      ##LOW, HIGH, RATE
      temp <- readBin(con, what="double", 3, size=4, endian="little")

      temp.LOW <- temp[1]
      temp.HIGH <- temp[2]
      temp.RATE <- temp[3]


      temp.TEMPERATURE<-readBin(con, what="integer", 1, size=2, endian="little")

      ##XCOORD, YCOORD, TOLDELAY, TOLON, TOLOFF
      temp <- readBin(con, what="integer", 5, size=2, endian="little")

      temp.XCOORD <- temp[1]
      temp.YCOORD <- temp[2]
      temp.TOLDELAY <- temp[3]
      temp.TOLON <- temp[4]
      temp.TOLOFF <- temp[5]


      ##POSITION
      temp.POSITION<-readBin(con, what="int", 1, size=1, endian="little")

      ##RUN
      temp.RUN<-readBin(con, what="int", 1, size=1, endian="little")

      ##TIME
      TIME_SIZE<-readBin(con, what="int", 1, size=1, endian="little")

      ##time size corrections for wrong time formats; set n to 6 for all values
      ##accoording the handbook of Geoff Duller, 2007
      TIME_SIZE<-6
      temp.TIME<-readChar(con, TIME_SIZE, useBytes=TRUE)


      ##DATE
      DATE_SIZE<-readBin(con, what="int", 1, size=1, endian="little")

      ##date size corrections for wrong date formats; set n to 6 for all values
      ##accoording the handbook of Geoff Duller, 2007
      DATE_SIZE<-6
      temp.DATE<-readChar(con, DATE_SIZE, useBytes=TRUE)


      ##SEQUENCE
      SEQUENCE_SIZE<-readBin(con, what="int", 1, size=1, endian="little")
      temp.SEQUENCE<-readChar(con, SEQUENCE_SIZE, useBytes=TRUE)

      #step forward in con
      if(8-SEQUENCE_SIZE>0){
        STEPPING<-readBin(con, what="raw", (8-c(SEQUENCE_SIZE)),size=1, endian="little")
      }


      ##USER
      USER_SIZE<-readBin(con, what="int", 1, size=1, endian="little")
      temp.USER<-readChar(con, USER_SIZE, useBytes=FALSE)

      #step forward in con
      if(8-c(USER_SIZE)>0){
        STEPPING<-readBin(con, what="raw", (8-c(USER_SIZE)), size=1, endian="little")
      }

      ##DTYPE
      temp.DTYPE<-readBin(con, what="int", 1, size=1, endian="little")

      ##IRR_TIME
      temp.IRR_TIME<-readBin(con, what="double", 1, size=4, endian="little")

      ##IRR_TYPE
      temp.IRR_TYPE<-readBin(con, what="int", 1, size=1, endian="little")

      ##IRR_UNIT
      temp.IRR_UNIT<-readBin(con, what="int", 1, size=1, endian="little")

      ##BL_TIME
      temp.BL_TIME<-readBin(con, what="double", 1, size=4, endian="little")

      ##BL_UNIT
      temp.BL_UNIT<-readBin(con, what="int", 1, size=1, endian="little")

      ##AN_TEMP, AN_TIME, NORM1, NORM2, NORM3, BG
      temp <- readBin(con, what="double", 6, size=4, endian="little")

      temp.AN_TEMP <- temp[1]
      temp.AN_TIME <- temp[2]
      temp.NORM1 <- temp[3]
      temp.NORM2 <- temp[4]
      temp.NORM3 <- temp[5]
      temp.BG <- temp[6]

      ##SHIFT
      temp.SHIFT<-readBin(con, what="integer", 1, size=2, endian="little")

      ##SAMPLE
      SAMPLE_SIZE<-readBin(con, what="int", 1, size=1, endian="little")
      temp.SAMPLE<-readChar(con, SAMPLE_SIZE, useBytes=TRUE) #however it should be set to 20

      #step forward in con
      if(20-c(SAMPLE_SIZE)>0){
        STEPPING<-readBin(con, what="raw", (20-c(SAMPLE_SIZE)), size=1, endian="little")
      }

      ##COMMENT
      COMMENT_SIZE<-readBin(con, what="int", 1, size=1, endian="little")
      temp.COMMENT<-readChar(con, COMMENT_SIZE, useBytes=TRUE) #set to 80 (manual)

      #step forward in con
      if(80-c(COMMENT_SIZE)>0){
        STEPPING<-readBin(con, what="raw", (80-c(COMMENT_SIZE)), size=1, endian="little")
      }

      ##LIGHTSOURCE, SET, TAG
      temp <- readBin(con, what="int", 3, size=1, endian="little")

      temp.LIGHTSOURCE <- temp[1]
      temp.SET <- temp[2]
      temp.TAG <- temp[3]


      ##GRAIN
      temp.GRAIN<-readBin(con, what="integer", 1, size=2, endian="little")

      ##LPOWER
      temp.LPOWER<-readBin(con, what="double", 1, size=4, endian="little")

      ##SYSTEMID
      temp.SYSTEMID<-readBin(con, what="integer", 1, size=2, endian="little")

      ##Unfortunately an inconsitent BIN-file structure forces a differenciation ...
      if(temp.VERSION == 03){

        ##RESERVED
        temp.RESERVED1<-readBin(con, what="raw", 36, size=1, endian="little")

        ##ONTIME, OFFTIME
        temp <- readBin(con, what="double", 2, size=4, endian="little")

        temp.ONTIME <- temp[1]
        temp.OFFTIME <- temp[2]


        ##Enable flags  #GateEnabled for v 06
        temp.ENABLE_FLAGS <- readBin(con, what="raw", 1, size=1, endian="little")
        temp.GATE_ENABLED <- temp.ENABLE_FLAGS

        ##ONGATEDELAY, OFFGATEDELAY
        temp <- readBin(con, what="double", 2, size=4, endian="little")

        temp.GATE_START <- temp[1]
        temp.GATE_STOP <- temp[2]

        ##RESERVED
        temp.RESERVED2<-readBin(con, what="raw", 1, size=1, endian="little")


      }else{

        ##RESERVED
        temp.RESERVED1<-readBin(con, what="raw", 20, size=1, endian="little")

        ##CURVENO
        temp.CURVENO <- readBin(con, what="integer", 1, size=2, endian="little")

        ##TIMETICK
        temp.TIMETICK <- readBin(con, what="double", 1, size=4, endian="little")

        ##ONTIME, STIMPERIOD
        temp <- readBin(con, what="integer", 2, size=4, endian="little")

        temp.ONTIME <- temp[1]
        temp.STIMPERIOD <- temp[2]

        ##GATE_ENABLED
        temp.GATE_ENABLED <- readBin(con, what="raw", 1, size=1, endian="little")

        ##ONGATEDELAY, OFFGATEDELAY
        temp <- readBin(con, what="double", 2, size=4, endian="little")

        temp.GATE_START <- temp[1]
        temp.GATE_END <- temp[2]
        temp.GATE_STOP <- temp.GATE_END

        ##PTENABLED
        temp.PTENABLED <- readBin(con, what="raw", 1, size=1, endian="little")

        ##RESERVED
        temp.RESERVED2<-readBin(con, what="raw", 10, size=1, endian="little")

      }

      #DPOINTS
      temp.DPOINTS<-readBin(con, what="integer", temp.NPOINTS, size=4, endian="little")


    }else{

      stop("[read_BIN2R()] Unsupported BIN/BINX-file version.")

    }

    #endif:format support
    ##END BIN FILE FORMAT SUPPORT
    ## ==========================================================================#

    #SET UNIQUE ID
    temp.ID <- temp.ID+1

    ##update progress bar
    if(txtProgressBar==TRUE){
      setTxtProgressBar(pb, seek(con,origin="current"))
    }
    ##set for equal values with different names
    if(exists("temp.GRAINNUMBER") == TRUE){temp.GRAIN <- temp.GRAINNUMBER}
    if(exists("temp.GRAIN") == TRUE){temp.GRAINNUMBER <- temp.GRAIN}

    if(exists("temp.LIGHTPOWER") == TRUE){temp.LPOWER <- temp.LIGHTPOWER}
    if(exists("temp.LPOWER") == TRUE){temp.LIGHTPOWER <- temp.LPOWER}

    temp.SEL <- if(temp.TAG == 1){TRUE}else{FALSE}

    ##replace values in the data.table with values
    results.METADATA[temp.ID, `:=` (
      ID = temp.ID,
      SEL = temp.SEL,
      VERSION = as.numeric(temp.VERSION),
      LENGTH = temp.LENGTH,
      PREVIOUS = temp.PREVIOUS,
      NPOINTS = temp.NPOINTS,
      RUN = temp.RUN,
      SET = temp.SET,
      POSITION = temp.POSITION,
      GRAIN = temp.GRAIN,
      GRAINNUMBER = temp.GRAINNUMBER,
      CURVENO = temp.CURVENO,
      XCOORD = temp.XCOORD,
      YCOORD = temp.YCOORD,
      SAMPLE = temp.SAMPLE,
      COMMENT = temp.COMMENT,
      SYSTEMID = temp.SYSTEMID,
      FNAME = temp.FNAME,
      USER = temp.USER,
      TIME = temp.TIME,
      DATE = temp.DATE,
      DTYPE = as.character(temp.DTYPE),
      BL_TIME = temp.BL_TIME,
      BL_UNIT = temp.BL_UNIT,
      NORM1 = temp.NORM1,
      NORM2 = temp.NORM2,
      NORM3 = temp.NORM3,
      BG = temp.BG,
      SHIFT = temp.SHIFT,
      TAG = temp.TAG,
      LTYPE = as.character(temp.LTYPE),
      LIGHTSOURCE = as.character(temp.LIGHTSOURCE),
      LPOWER = temp.LPOWER,
      LIGHTPOWER = temp.LIGHTPOWER,
      LOW = temp.LOW,
      HIGH = temp.HIGH,
      RATE = temp.RATE,
      TEMPERATURE = temp.TEMPERATURE,
      MEASTEMP = temp.MEASTEMP,
      AN_TEMP = temp.AN_TEMP,
      AN_TIME = temp.AN_TIME,
      TOLDELAY = temp.TOLDELAY,
      TOLON = temp.TOLON,
      TOLOFF = temp.TOLOFF,
      IRR_TIME = temp.IRR_TIME,
      IRR_TYPE = temp.IRR_TYPE,
      IRR_UNIT = temp.IRR_UNIT,
      IRR_DOSERATE = temp.IRR_DOSERATE,
      IRR_DOSERATEERR = temp.IRR_DOSERATEERR,
      TIMESINCEIRR = temp.TIMESINCEIRR,
      TIMETICK = temp.TIMETICK,
      ONTIME = temp.ONTIME,
      OFFTIME = temp.OFFTIME,
      STIMPERIOD = temp.STIMPERIOD,
      GATE_ENABLED = as.numeric(temp.GATE_ENABLED),
      ENABLE_FLAGS = as.numeric(temp.ENABLE_FLAGS),
      GATE_START = temp.GATE_START,
      GATE_STOP = temp.GATE_STOP,
      PTENABLED = as.numeric(temp.PTENABLED),
      DTENABLED = as.numeric(temp.DTENABLED),
      DEADTIME = temp.DEADTIME,
      MAXLPOWER = temp.MAXLPOWER,
      XRF_ACQTIME = temp.XRF_ACQTIME,
      XRF_HV = temp.XRF_HV,
      XRF_CURR = temp.XRF_CURR,
      XRF_DEADTIMEF = temp.XRF_DEADTIMEF,
      DETECTOR_ID = temp.DETECTOR_ID,
      LOWERFILTER_ID = temp.LOWERFILTER_ID,
      UPPERFILTER_ID = temp.UPPERFILTER_ID,
      ENOISEFACTOR = temp.ENOISEFACTOR,
      SEQUENCE = temp.SEQUENCE


    )]

    results.DATA[[temp.ID]] <- temp.DPOINTS

    results.RESERVED[[temp.ID]][[1]] <- temp.RESERVED1
    results.RESERVED[[temp.ID]][[2]] <- temp.RESERVED2

    ##BREAK
    ##stop loop if record limit is reached
    if (!is.null(n.records)) {
      if (n.records == temp.ID) {
        break()
      }

    }

    ##remove some unwanted objects
    rm(temp.GRAINNUMBER)
    rm(temp.GRAIN)


  }#endwhile::end lopp


  ##close con
  close(con)

  ##close
  if(txtProgressBar==TRUE){close(pb)}

  ##output
  cat(paste("\t >> ",temp.ID," records have been read successfully!\n\n", sep=""))

  # Further limitation --------------------------------------------------------------------------
  if(!is.null(position)){

    ##check whether the position is valid at all
    if (all(position %in% results.METADATA$POSITION)) {

      results.METADATA <- results.METADATA[which(results.METADATA$POSITION %in% position),]
      results.DATA <- results.DATA[results.METADATA$ID]

        ##re-calculate ID ... otherwise it will not match
        results.METADATA$ID <- 1:length(results.DATA )

        ##show a message
        message("[read_BIN2R()] The record index has been recalculated!")


    }else{
      valid.position <-
        paste(unique(results.METADATA$POSITION), collapse = ", ")
      warning(
        paste0(
          "Position limitation omitted. At least one position number is not valid, valid position numbers are: ", valid.position
        )
      )
    }

  }

  ##check for duplicated entries and remove them if wanted
  duplication.check <- suppressWarnings(which(c(0, vapply(
    2:length(results.DATA),
    FUN = function(x) {
      all(results.DATA[[x - 1]] == results.DATA[[x]])
    },
    FUN.VALUE = 1
  )) == 1))

  if(length(duplication.check) != 0){

    if(duplicated.rm){

      ##remove records
      results.METADATA <- results.METADATA[-duplication.check,]
      results.DATA[duplication.check] <- NULL

      ##message
      message(paste0("[read_BIN2R()] duplicated record(s) detected and removed: ", paste(duplication.check, collapse = ", ")))

    }else{

      warning(paste0("[read_BIN2R()] duplicated record(s) detected: ",
                     paste(duplication.check, collapse = ", "), " You should consider 'duplicated.rm = TRUE'." ))

    }

  }

  ##produce S4 object for output
  object <- set_Risoe.BINfileData(
    METADATA = results.METADATA,
    DATA = results.DATA,
    .RESERVED =  results.RESERVED)

  # Convert Translation Matrix Values ---------------------------------------

  if (!show.raw.values) {
    ##LIGHTSOURCE CONVERSION
    object@METADATA[, "LIGHTSOURCE"] <-
      unname(LIGHTSOURCE.lookup[object@METADATA[, "LIGHTSOURCE"]])

    ##LTYPE CONVERSION
    object@METADATA[, "LTYPE"] <-
      unname(LTYPE.lookup[object@METADATA[, "LTYPE"]])

    ##DTYPE CONVERSION
    object@METADATA[, "DTYPE"] <-
      unname(DTYPE.lookup[object@METADATA[, "DTYPE"]])

        ##CHECK for oddly set LTYPES, this may happen in old BIN-file versions
        if (object@METADATA[1, "VERSION"] == 3) {
          object@METADATA[, "LTYPE"] <-
            sapply(1:length(object@METADATA[, "LTYPE"]), function(x) {
              if (object@METADATA[x, "LTYPE"] == "OSL" &
                  object@METADATA[x, "LIGHTSOURCE"] == "IR diodes/IR Laser") {
                return("IRSL")

              } else{
                return(object@METADATA[x, "LTYPE"])

              }

            })

        }

    ##TIME CONVERSION, do not do for odd time formats as this could cause problems during export
    if (TIME_SIZE == 6) {
      object@METADATA[, "TIME"] <-
        sapply(1:length(object@METADATA[, "TIME"]), function(x) {
          format(strptime(as.character(object@METADATA[x, "TIME"]), "%H%M%S"), "%H:%M:%S")

        })
    }

  }


  # Fast Forward --------------------------------------------------------------------------------
  ## set fastForward to TRUE if one of this arguments is used
  if(any(names(list(...)) %in% names(formals(Risoe.BINfileData2RLum.Analysis))[-1]) &
     fastForward == FALSE) {
    fastForward <- TRUE
    warning("[read_BIN2R()] automatically reset 'fastForward = TRUE'")

  }

  ##return values
  ##with fast fastForward they will be converted directly to a list of RLum.Analysis objects
  if(fastForward){
     object <- Risoe.BINfileData2RLum.Analysis(object, ...)


     ##because we expect a list
     if(!is(object, "list")){
       object <- list(object)

     }


  }


   return(object)


}

## ---- DEPRECATED GENERICS
# .Deprecated in package version 0.5.0
# .Defunct in 0.5.1
# Removed in 0.6.0
#' @noRd
#' @export
readBIN2R <- function(...) {
  .Defunct("read_BIN2R")
  read_BIN2R(...)
}

