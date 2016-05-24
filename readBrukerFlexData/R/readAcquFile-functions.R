## Copyright 2010-2014 Sebastian Gibb
## <mail@sebastiangibb.de>
##
## This file is part of readBrukerFlexData for R and related languages.
##
## readBrukerFlexData is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## readBrukerFlexData is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with readBrukerFlexData. If not, see <http://www.gnu.org/licenses/>

#' Reads an acqu file.
#'
#' This function reads constants, calibrations values and a lot of more from
#' \emph{acqu} files.
#'
#' @param file \code{character}, path to corresponding fid file
#'  (e.g. \code{Pankreas_HB_L_061019_A10/0_a19/1/1SLin/fid})
#' @param verbose \code{logical}, print verbose messages?
#'
#' @return
#'  a \code{list} with metadata
#'
#' @seealso \code{\link[readBrukerFlexData]{readBrukerFlexFile}},
#' @rdname readAcquFile
#' @keywords internal IO
#'
##  We have to import the following data to calculating mass:
##  $BYTORDA: endianness, 0==little, 1==big
##      => metaData$byteOrder
##  $TD: total number of measured time periods
##      => metaData$number
##  $DELAY: first measured intensity after $DELAY ns
##      => metaData$timeDelay
##  $DW: ns between measured time periods
##      => metaData$timeDelta
##  $ML1: calibration constant
##      => metaData$calibrationConstants[1]
##  $ML2: calibration constant
##      => metaData$calibrationConstants[2]
##  $ML3: calibration constant
##      => metaData$calibrationConstants[3]
##
##  if we want to use High Precision Calibration (HPC), we need:
##  $HPClBHi: upper mass threshold
##      => metaData$hpcLimits["maxMass"]
##  $HPClBLo: lower mass threshold
##      => metaData$hpcLimits["minMass"]
##  $HPClOrd: polynomial order
##      => metaData$hpcOrder
##  $HPClUse: maybe using of HPC? (seems always be "yes" in our test data)
##      => metaData$hpcUse
##  $HPCStr: polynomial coefficients in a string
##      => metaData$hpcCoefficients
##
##  we try to import [optional]:
##  DATATYPE
##  SPECTROMETER/DATASYSTEM
##      => metaData$dataSystem
##  .SPECTROMETER TYPE
##      => metaData$spectrometerType
##  .INLET
##      => metaData$inlet
##  .IONIZATION MODE
##      => metaData$ionizationMode
##  $DATE
##      => metaData$date
##  $ACQMETH
##      => metaData$acquisitionMethod
##  $AQ_DATE
##      => metaData$acquisitionDate
##  $AQ_mod
##      => metaData$acquisitionMode
##  $AQOP_mod
##      => metaData$acquisitionOperatorMode
##      => metaData$tofMode (replaces file path based method)
##  $ATTEN
##      => metaData$laserAttenuation
##  $COM[1:4]
##      => metaData$comments
##  $DEFLON
##      => metaData$deflection
##  $DIGTYP
##      => metaData$digitizerType
##  $DPCAL1
##      => metaData$deflectionPulserCal1
##  $DPMASS
##      => metaData$deflectionPulserMass
##  $FCVer
##      => metaData$flexControlVersion
##  $ID_raw
##      => metaData$id
##  $INSTRUM
##      => metaData$instrument
##  $InstrID
##      => metaData$instrumentId
##  $InstrTyp
##      => metaData$instrumentType
##  $Masserr
##      => metaData$massError
##  $NoSHOTS: number of laser shots
##      => metaData$laserShots
##  $SPType
##      => metaData$spectrumType
##  $PATCHNO: sample position on target
##      => metaData$patch
##  $PATH: original file path (on Bruker flex series controller PC)
##      => metaData$path
##  $REPHZ
##      => metaData$laserRepetition
##  $SPOTNO: same as $PATCHNO (in older files often empty, that's why we use
##      $PATCHNO instead)
##      => metaData$spot
##  $TgIDS: target ids
##      => metaData$targetIdString
##  $TgCount: number of measurement with this target
##      => metaData$targetCount
##  $TgSer: target serial number
##      => metaData$targetSerialNumber
##  $TgTyp: target type number
##      => metaData$targetTypeNumber
##
##  import from file path:
##  full current path to fid file:
##      => metaData$fidFile
##  sample name
##      => metaData$sampleName
##
.readAcquFile <- function(fidFile, verbose=FALSE) {
  acquFile <- sub(pattern="fid$", x=fidFile, replacement="acqu")

  if (verbose) {
    message("Reading metadata from ", sQuote(acquFile), " ...")
  }

  if (!file.exists(acquFile)) {
    stop("File ", sQuote(acquFile), " doesn't exists!")
  }

  con <- file(acquFile, "rt")
  acquLines <- readLines(con, n=-1)
  close(con)

  ## collect data
  metaData <- list()

  ## endianness
  isBigEndian <- as.integer(.grepAcquValue("##\\$BYTORDA=", acquLines)) == 1
  metaData$byteOrder <- ifelse(isBigEndian, "big", "little")

  ## obligate
  metaData$number <- as.double(.grepAcquValue("##\\$TD=", acquLines))
  metaData$timeDelay <- .grepAcquDoubleValue("##\\$DELAY=", acquLines)
  metaData$timeDelta <- .grepAcquDoubleValue("##\\$DW=", acquLines)
  metaData$calibrationConstants <-
    c(c1=.grepAcquDoubleValue("##\\$ML1=", acquLines),
      c2=.grepAcquDoubleValue("##\\$ML2=", acquLines),
      c3=.grepAcquDoubleValue("##\\$ML3=", acquLines))

  ## obligate HPC
  metaData$hpcLimits <-
    c(minMass=.grepAcquDoubleValue("##\\$HPClBLo=", acquLines),
      maxMass=.grepAcquDoubleValue("##\\$HPClBHi=", acquLines))
  metaData$hpcOrder <- as.double(.grepAcquValue("##\\$HPClOrd=", acquLines))
  metaData$hpcUse <-
    as.logical(.grepAcquValue("##\\$HPClUse=", acquLines) == "yes")

  ## was HPC involved?  metaData$hpcUse seems to be always true
  isHPCused <- isTRUE(metaData$hpcUse &&
                      metaData$hpcLimits["maxMass"] > 0 &&
                      metaData$hpcLimits["minMass"] > 0 &&
                      metaData$hpcOrder > 0)

  if (isHPCused) {
    hpcStr <- .grepAcquValue("##\\$HPCStr=", acquLines)
    hpcConstants <- .extractHPCConstants(hpcStr)
    metaData$hpcCoefficients <- hpcConstants$coefficients
    metaData$hpcCalibrationConstant0 <- hpcConstants$calibrationConstant0
    metaData$hpcCalibrationConstant2 <- hpcConstants$calibrationConstant2
  }

  ## obligate LIFT
  metaData$lift <- c(.grepAcquDoubleValue("##\\$Lift1=", acquLines),
                     .grepAcquDoubleValue("##\\$Lift2=", acquLines))
  metaData$tlift <- .grepAcquDoubleValue("##\\$TLift=", acquLines)

  ## optional
  metaData$dataType <- .grepAcquValue("##DATATYPE=", acquLines)
  metaData$dataSystem <- .grepAcquValue("##SPECTROMETER/DATASYSTEM=", acquLines)
  metaData$spectrometerType <-
    .grepAcquValue("##.SPECTROMETER TYPE=", acquLines)
  metaData$inlet <- .grepAcquValue("##.INLET=", acquLines)
  metaData$ionizationMode <- .grepAcquValue("##.IONIZATION MODE=", acquLines)
  metaData$date <- .grepAcquValue("##\\$DATE=", acquLines)


  metaData$acquisitionMethod <- .grepAcquValue("##\\$ACQMETH=", acquLines)
  metaData$acquisitionDate <- .grepAcquValue("##\\$AQ_DATE=", acquLines)
  aq_mod <- .grepAcquValue("##\\$AQ_mod=", acquLines)
  if (length(aq_mod)) {
    metaData$acquisitionMode <- switch(aq_mod,
      "0" = { "qf" },
      "1" = { "qsim" },
      "2" = { "qseq" },
      { aq_mod }
    )
  }

  aqop <- .grepAcquValue("##\\$AQOP_m=", acquLines)
  if (length(aqop)) {
    metaData$tofMode  <- switch(aqop,
      "0" = { "LINEAR" },
      "1" = { "REFLECTOR" },
      { aqop }
    )
  }

  metaData$acquisitionOperatorMode <- metaData$tofMode

  metaData$laserAttenuation <- .grepAcquDoubleValue("##\\$ATTEN=", acquLines)

  metaData$comments <- .grepAcquValue("##\\$CMT.*=", acquLines)

  metaData$deflection <-
    as.logical(.grepAcquValue("##\\$DEFLON=", acquLines) == "yes")

  digtyp  <- .grepAcquValue("##\\$DIGTYP=", acquLines)
  if (length(digtyp)) {
    metaData$digitizerType <- switch(digtyp,
        "0" = { "unknown" },
        "1" = { "Lecroy LSA1000" },
        "2" = { "Acqiris DP105" },
        "3" = { "Acqiris DP110" },
        "4" = { "Acqiris DP211" },
        "5" = { "Acqiris DP240" },
        "6" = { "Acqiris AP200" },
        "7" = { "Acqiris AP240" },
        "8" = { "Acqiris DC440" },
        "9" = { "Acqiris DC282" },
       "10" = { "Acqiris Unknown subtype" },
       "11" = { "Gage" },
       "12" = { "Simulator" },
       "13" = { "Lecroy WaveRunner" },
       "14" = { "Acqiris U1084A" },
       "15" = { "NI 5154" },
       "16" = { "LeCroy LSA2000" },
       "17" = { "Acqiris DP1400" },
       "18" = { "NI 5155" },
       "19" = { "Bruker BD0G5" },
        { digtyp }
    )
  }

  metaData$deflectionPulserCal1 <-
    .grepAcquDoubleValue("##\\$DPCAL1=", acquLines)
  metaData$deflectionPulserMass <-
    .grepAcquDoubleValue("##\\$DPMASS=", acquLines)
  metaData$flexControlVersion <- .grepAcquValue("##\\$FCVer=", acquLines)
  metaData$id <- .grepAcquValue("##\\$ID_raw=", acquLines)

  metaData$instrument <- .grepAcquValue("##\\$INSTRUM=", acquLines)
  metaData$instrumentId <- .grepAcquValue("##\\$InstrID=", acquLines)
  metaData$instrumentType <- .grepAcquValue("##\\$InstTyp=", acquLines)

  metaData$massError <- .grepAcquDoubleValue("##\\$Masserr=", acquLines)

  metaData$laserShots <- as.double(.grepAcquValue("##\\$NoSHOTS=", acquLines))

  if (metaData$laserShots == 0) {
    warning("File ", sQuote(fidFile), " seems to be empty because ",
            "no laser shots applied to this sample.")
  }

  metaData$patch <- .grepAcquValue("##\\$PATCHNO=", acquLines)

  ## imaging data
  if (length(metaData$patch) &&
      grepl(pattern="(R[0-9]+)?X[0-9]+Y[0-9]+", x=metaData$patch,
            ignore.case=TRUE)) {
    rx <- gregexpr(pattern="[XY][0-9]+", text=metaData$patch)[[1]]
    pos <- substring(metaData$patch, rx+1, rx+attr(rx, "match.length")-1)

    if (length(pos) == 2) {
      pos <- as.double(pos)
      metaData$imaging <- list(pos=c(x=pos[1], y=pos[2]))
    }
  }

  metaData$path <- .grepAcquValue("##\\$PATH=", acquLines)
  metaData$laserRepetition <- .grepAcquDoubleValue("##\\$REPHZ=", acquLines)
  metaData$spot <- .grepAcquValue("##\\$SPOTNO=", acquLines)

  sptype <- .grepAcquValue("##\\$SPType=", acquLines)
  if (length(sptype)) {
    metaData$spectrumType <- switch(sptype,
      "0" = { "TOF" },
      "1" = { "PSD" },
      "2" = { "LIFT" },
      "3" = { "PSDSegment" },
      { sptype }
    )
  }

  metaData$targetCount <- as.double(.grepAcquValue("##\\$TgCount", acquLines))
  metaData$targetIdString <- .grepAcquValue("##\\$TgIDS", acquLines)
  metaData$targetSerialNumber <- .grepAcquValue("##\\$TgSer", acquLines)
  metaData$targetTypeNumber <- .grepAcquValue("##\\$TgTyp", acquLines)

  metaData$file <- fidFile

  metaData$sampleName <- .sampleName(fidFile)
  metaData$fullName <- paste(metaData$sampleName, metaData$patch, sep=".")
  metaData$name <- metaData$fullName

  return(metaData)
}

