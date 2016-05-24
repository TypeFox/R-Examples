## Copyright 2010-2013 Sebastian Gibb
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

#' Reads mass spectrometry data in Bruker Daltonics XMASS format.
#'
#' This function reads mass spectrometry data in Bruker Daltonics XMASS format
#' used by Bruker Daltonics mass spectrometer of *flex series (autoflex,
#' microflex, ultraflex).
#'
#' @details
#' \code{readBrukerFlexFile} has to import the following data to calculating
#' mass from \emph{acqu} file:\cr
#'
#' \tabular{lll}{
#' acqu-value \tab becomes metaData \tab description \cr
#' $BYTORDA \tab metaData$byteOrder \tab endianness of fid file \cr
#' $TD \tab metaData$number \tab total number of measured time periods \cr
#' $DELAY \tab metaData$timeDelay \tab first measured intensity after
#'   \emph{metaData$timeDelay} ns \cr
#' $DW \tab metaData$timeDelta \tab ns between measured time periods \cr
#' $ML1 \tab metaData$calibrationConstants[1] \tab mass calibration constant \cr
#' $ML2 \tab metaData$calibrationConstants[2] \tab mass calibration constant \cr
#' $ML3 \tab metaData$calibrationConstants[3] \tab mass calibration constant \cr
#' }
#'
#' If High Precision Calibration (HPC) is used, \code{readBrukerFlexFile} needs:
#' \tabular{lll}{
#' acqu-value \tab becomes metaData \tab description \cr
#' $HPClBHi \tab metaData$hpc$limits[\dQuote{maxMass}] \tab upper mass
#'  threshold \cr
#' $HPClBLo \tab metaData$hpc$limits[\dQuote{minMass}] \tab lower mass
#'  threshold \cr
#' $HPClOrd \tab metaData$hpc$order \tab polynomial order \cr
#' $HPClUse \tab metaData$hpc$use \tab maybe using of HPC? (seems to be always
#' \dQuote{yes} in our test data) \cr
#' $HPCStr \tab metaData$hpc$coefficients \tab polynomial coefficients in a
#' string \cr
#' }
#'
#' \code{readBrukerFlexFile} tries also to import [optional]: \cr
#' \tabular{lll}{
#' acqu-value \tab becomes metaData \tab description \cr
#' DATATYPE \tab metaData$dataType \tab e.g CONTINUOUS MASS SPECTRUM \cr
#' SPECTROMETER/DATASYSTEM \tab metaData$dataSystem \tab
#'  e.g. Bruker Flex Series \cr
#' .SPECTROMETER TYPE \tab metaData$spectrometerType \tab e.g. TOF \cr
#' .INLET \tab metaData$inlet \tab DIRECT \cr
#' .IONIZATION MODE \tab metaData$ionizationMode \tab e.g. LD+ \cr
#' $DATE \tab metaData$date \tab same as $AQ_DATE but often only "0" \cr
#' $ACQMETH \tab metaData$acquisitionMethod \tab path to method file \cr
#' $AQ_DATE \tab metaData$acquisitionDate \tab acquisition date \cr
#' $AQ_mod \tab metaData$acquisitionMode \tab acquisition mode \cr
#' $AQOP_m \tab metaData$acquisitionOperatorMode, metaData$tofMode \tab
#'  LINEAR / REFLECTOR \cr
#' $ATTEN \tab metaData$laserAttenuation \tab laser beam attenuation \cr
#' $CMT[1:4] \tab metaData$comments \tab comments \cr
#' $DEFLON \tab metaData$deflection \tab deflection ON/OFF \cr
#' $DIGTYP \tab metaData$digitizerType \tab type of digitizer \cr
#' $DPCAL1 \tab metaData$deflectionPulserCal1 \tab deflection pulser cal 1 \cr
#' $DPMASS \tab metaData$deflectionPulserMass \tab deflection pulser mass \cr
#' $FCVer \tab metaData$flexControlVersion \tab Version of
#'  Bruker Daltonics FlexControl software \cr
#' $ID_raw \tab metaData$id \tab spectrum id \cr
#' $INSTRUM \tab metaData$instrument \tab e.g. AUTOFLEX \cr
#' $InstrID \tab metaData$instrumentId \tab ID of mass spectrometer \cr
#' $InstTyp \tab metaData$instrumentType \tab instrument type \cr
#' $Lift1 \tab metaData$lift[1] \tab LIFT constant? \cr
#' $Lift2 \tab metaData$lift[2] \tab LIFT constant? \cr
#' $Masserr \tab metaData$massError \tab initial mass error in ppm \cr
#' $NoSHOTS \tab metaData$laserShots \tab number of applied laser shots \cr
#' $PATCHNO \tab metaData$patch \tab sample postion on target \cr
#' $PATH \tab metaData$path \tab original file path
#'  (on Bruker *flex series controller PC) \cr
#' $REPHZ \tab metaData$laserRepetition \tab laser repetition rate in Hz \cr
#' $SPOTNO \tab metaData$spot \tab same as \emph{$PATCHNO}
#'  (in older files often empty, that's why \code{readBrukerFlexFile}
#'  uses \emph{$PATHNO} instead) \cr
#' $SPType \tab metaData$spectrumType \tab e.g. TOF \cr
#' $TgIDS \tab metaData$target$id \tab target ids \cr
#' $TgCount \tab metaData$target$count \tab number of measurements
#'  with this target \cr
#' $TgSer \tab metaData$target$serialNumber \tab target serial number \cr
#' $TgTyp \tab metaData$target$typeNumber \tab target type number \cr
#' $TLift \tab metaData$tlift \tab LIFT constant? \cr
#' }
#'
#' import from file path:
#' \tabular{lll}{
#' value \tab becomes metaData \tab description \cr
#' full current path to fid file \tab metaData$file \tab
#'  path on local machine \cr
#' sample name \tab metaData$sampleName \tab - \cr
#' }
#'
#' \code{filterZeroIntensities}: Change default value is \bold{not recommended}!
#' If \code{TRUE} all intensities equal zero are removed.
#' This parameter exists only to be compatible to
#' Bruker Daltonics CompassXport's mzXML export function. For details see:
#' \sQuote{Release Notes for CompassXport 3.0.3},
#' cap. 6 \sQuote{Filtering of Zero Intensities}:
#' \dQuote{Bruker Daltonics' Acquisition Software will compress Analysis raw
#' data. To save on operation time and to keep export file sizes small,
#' CompassXport 3.0.3 will filter out zero (0.0) intensities
#' when exporting to mzXML or mzData \ldots}
#'
#' \code{keepNegativeIntensities}: Change default value is
#' \bold{not recommended}! If \code{TRUE} negative intensity values are not
#' replaced by zero. This parameter exists only to be compatible to
#' Bruker Daltonics CompassXport.
#'
#'
#' @param fidFile \code{character}, path to \emph{fid} file which should be
#'  read.
#' @param removeMetaData, \code{logical}, to calculate mass data a lot of meta
#'  data are needed. To save memory they could be deleted after calculation.
#' @param useHpc \code{logical}, should Bruker Daltonics' High Precision
#'  Calibration be used if available?
#'  (see also: \code{\link[readBrukerFlexData]{.hpc}})
#' @param filterZeroIntensities \code{logical}, don't change it.
#'  If \code{TRUE} all intensities equal \code{0.0} are removed.
#'  (see also: \sQuote{Details} section)
#' @param keepNegativeIntensities \code{logical}, don't change it.
#'  If \code{FALSE} all intensities less than zero are replaced by
#'  zero.
#'  (see also: \sQuote{Details} section)
#' @param verbose \code{logical}, print verbose messages?
#'
#' @return
#'  A \code{list} of spectra and metadata.
#'
#'  \itemize{
#'     \item{\code{spectrum$mass}: }{A vector of calculated mass.}
#'     \item{\code{spectrum$tof}: }{A vector of time-of-flight data.}
#'     \item{\code{spectrum$intensity}: }{A vector of intensity values.}
#'     \item{\code{metaData}: }{A list of metaData depending on read spectrum.}
#' }
#' @export
#' @seealso
#'  \url{https://github.com/sgibb/readBrukerFlexData/wiki},
#'  \code{\link[MALDIquantForeign]{importBrukerFlex}},
#'  \code{\link[readBrukerFlexData]{readBrukerFlexDir}},
#'  \code{\link[readBrukerFlexData]{.hpc}}
#' @keywords IO
#' @rdname readBrukerFlexFile
#' @examples
#' ## load library
#' library("readBrukerFlexData")
#'
#' ## get examples directory
#' exampleDirectory <- system.file("Examples", package="readBrukerFlexData")
#'
#' ## read example spectrum
#' spec <- readBrukerFlexFile(file.path(exampleDirectory,
#'   "2010_05_19_Gibb_C8_A1/0_A1/1/1SLin/fid"))
#'
#' ## print metaData
#' print(spec$metaData)
#'
#' ## plot spectrum
#' plot(spec$spectrum$mass, spec$spectrum$intensity, type="l", col="red")
#'
readBrukerFlexFile <- function(fidFile, removeMetaData=FALSE, useHpc=TRUE,
                               filterZeroIntensities=FALSE,
                               keepNegativeIntensities=FALSE,
                               verbose=FALSE) {
  if (verbose) {
    message("Reading spectrum from ", sQuote(fidFile), " ...")
  }

  if (!file.exists(fidFile)) {
    stop("File ", sQuote(fidFile), " doesn't exists!")
  }

  if (file.info(fidFile)$isdir) {
    stop("Not a fid file! ", sQuote(fidFile), " is a directory.")
  }

  ## try to get absolute file path
  fidFile <- normalizePath(fidFile)

  ## read metadata
  metaData <- .readAcquFile(fidFile=fidFile, verbose=verbose)

  ## read peak intensities
  intensity <- .readFidFile(fidFile, metaData$number, metaData$byteOrder)

  ## calculate tof of metadata
  tof <- as.double(metaData$timeDelay +
                   ((0:(metaData$number-1)) * metaData$timeDelta))

  ## replace negative intensity values by zero
  isNegative <- which(intensity < 0)
  if (!keepNegativeIntensities && length(isNegative)) {
    intensity[isNegative] <- 0
  }

  ## remove times which have no intensity (intensity == 0), e.g. generated by
  ## REFLECTOR mode
  ## for details see "Release Notes for CompassXport 3.0.3"
  ## cap. 6 "Filtering of Zero Intensities"
  ## "Bruker Daltonics' Acquisition Software will compress Analysis raw
  ## data. To save on operation time and to keep export file sizes small,
  ## CompassXport 3.0.3 will filter out zero (0.0) intensities
  ## when exporting to mzXML or mzData ..."
  if (filterZeroIntensities) {
    notNull <- which(intensity > 0)
    intensity <- intensity[notNull]
    tof <- tof[notNull]
  }

  ## calculate mass of TOFs
  mass <- .tof2mass(tof,
                    metaData$calibrationConstants[[1]],
                    metaData$calibrationConstants[[2]],
                    metaData$calibrationConstants[[3]])

  ## TODO: fix equations in .hpc

  ## was HPC involved?
  ## metaData$hpcUse seems to be always true
  isHPCused <- isTRUE(useHpc & length(metaData$hpcCoefficients))

  if (isHPCused) {
    ## TODO: fix equations in .hpc and remove the following warning
    warning("The spectrum file ", sQuote(fidFile), " uses HPC. ",
            "HPC isn't fully supported by readBrukerFlexFile. ",
            "Please see ", dQuote("?.hpc"), " for details.")
    mass <- .hpc(mass=mass,
                 minMass=metaData$hpcLimits["minMass"],
                 maxMass=metaData$hpcLimits["maxMass"],
                 hpcCoefficients=metaData$hpcCoefficients)
  }

  ## TODO: add LIFT support

  ## was LIFT involved?
  isLIFTused <- isTRUE(all(metaData$lift != 0) || metaData$tlift != 0 ||
                       metaData$spectrumType == "LIFT")

  if (isLIFTused) {
    ## TODO: add LIFT support
    warning("The spectrum file ", sQuote(fidFile), " uses LIFT.\n",
            "LIFT isn't supported by readBrukerFlexFile. ",
            "Could not convert time-of-flight values into mass!")
    mass <- tof
  }

  spectrum <- list(tof=tof, mass=mass, intensity=intensity)

  if (!removeMetaData) {
    return(list(spectrum=spectrum, metaData=metaData))
  } else {
    return(list(spectrum=spectrum, metaData=list(file=metaData$file)))
  }
}

