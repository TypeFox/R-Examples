## Copyright 2010-2012 Sebastian Gibb
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

#' Extracts High Precision Calibration coefficients.
#'
#' This is a helper function to extract coefficients and constants from
#' metaData$hpcStr.
#'
#' @param hpcStr metaData$hpcStr, which store coefficents
#' \preformatted{
#' hpcStr <- "V1.0CHPCData  Order 10 vCoeff V1.0VectorDouble 11
#'            -0.48579224953906053 0.0009361303203700988
#'            -6.92711401708155e-008 -1.0992953299897006e-009
#'            1.1718229914003113e-012 -5.392578762547374e-016
#'            9.0176664604755316e-020 1.9704001597871883e-023
#'            -1.1794161284667635e-026 2.0351573912658823e-030
#'            -1.2617853301428769e-034
#'            c2 -0.046701600316874939
#'            c0 237.64781433281422
#'            minMass 736.50266799999997
#'            maxMass 3698.6377320000001 bUse 1 endCHPCData"
#' }
#'
#' @return
#'  a \code{list}:
#'  \itemize{
#'    \item{\code{hpcConstants$coefficients}: }{double, vector of coefficents}
#'    \item{\code{hpcConstants$calibrationContant0}: }{c0}
#'    \item{\code{hpcConstants$calibrationContant2}: }{c2}
#'  }
#'
#' @seealso \code{\link[readBrukerFlexData]{.hpc}}
#' @rdname extractHPCConstants
#' @keywords internal
#'
.extractHPCConstants <- function(hpcStr) {
  tmpLine <- strsplit(x=hpcStr, split=" ")[[1]]
  ## remove emtpy elements
  tmpLine <- tmpLine[tmpLine != ""]

  hpcConstants <- list()

  ## extract only coefficients
  hpcConstants$coefficients <-
    as.double(tmpLine[(which(tmpLine == "V1.0VectorDouble") + 2):(which(tmpLine == "c2") - 1)])

  hpcConstants$calibrationConstant2 <-
    as.double(tmpLine[which(tmpLine == "c2") + 1])
  hpcConstants$calibrationConstant0 <-
    as.double(tmpLine[which(tmpLine == "c0") + 1])

  return(hpcConstants)
}

#' High Precision Calibration
#'
#' Only basic support (not 100\% identical results) for Bruker Daltonics' HPC.
#' HPC stands for \bold{H}igh \bold{P}recision \bold{C}alibration.\cr
#' This is an internal function and should normally not used by the user.
#'
#' @param mass \code{double}, mass calculated traditionally.
#' @param minMass \code{double}, lower Threshold for HPC. HPC is only defined
#'  for a range of mass.
#' @param maxMass \code{double}, upper Threshold for HPC. HPC is only defined
#'  for a range of mass.
#' @param hpcCoefficients \code{doubles}, coefficients needed by the HPC
#'  algorithm (see also:
#'  \code{\link[readBrukerFlexData]{.extractHPCConstants}})
#' @return A vector of HPC corrected mass (\code{double}).
#'
#' @details
#' Bruker Daltonics doesn't explain how HPC works. All formula are results of
#' \dQuote{trial and error}. That is why mass calculated by \code{.hpc}
#' differs little from original HPC mass.
#' (In example file 214 of 24860 mass are incorrect;
#' deviations: min: 6.103515625e-05, max: 0.02935791015625.) \cr
#' In the manual of mass spectrometry instruments of Bruker Daltonics machines
#' the *flex series you can find an article about HPC principles: \cr
#'  Gobom, J. and Mueller, M. and Egelhofer V. and Theiss, D. and
#'  Lehrach, H. and Nordhoff, E. (2002) \cr
#'  \dQuote{A Calibration Method That Simplifies and Improves Accurate
#'   Determination of Peptide Molecular mass by MALDI-TOF MS.},
#'  \emph{Anal Chem} \bold{74}: 3915-3923 \cr
#'  \url{http://www.ncbi.nlm.nih.gov/pubmed/12175185}
#'
#' @note
#' Please note that .hpc is not correct! You have been warned.
#'
#' @references
#'  Gobom, J. and Mueller, M. and Egelhofer V. and Theiss, D. and
#'  Lehrach, H. and Nordhoff, E. (2002) \cr
#'  \dQuote{A Calibration Method That Simplifies and Improves Accurate
#'   Determination of Peptide Molecular mass by MALDI-TOF MS.},
#'  \emph{Anal Chem} \bold{74}: 3915-3923 \cr
#'  \url{http://www.ncbi.nlm.nih.gov/pubmed/12175185}
#'
#' @seealso
#'     \code{\link[readBrukerFlexData]{readBrukerFlexDir}},
#'     \code{\link[readBrukerFlexData]{readBrukerFlexFile}},
#'     \code{\link[readBrukerFlexData]{.double2singlePrecision}}
#'
#' @keywords internal
#' @rdname hpc
#' @examples
#' ## load library
#' library("readBrukerFlexData")
#'
#' ## get examples directory
#' exampleDirectory <- system.file("Examples", package="readBrukerFlexData")
#'
#' ## read example spectra
#' ## please note: filterZeroIntensities=TRUE is used for compatibility reasons.
#' ##              You should NOT use it!
#' noHpcSpec <- readBrukerFlexFile(file.path(exampleDirectory,
#'     "hpc/fid/0_A20/1/1SRef/fid"), filterZeroIntensities=TRUE, useHpc=FALSE)
#' hpcSpec <- readBrukerFlexFile(file.path(exampleDirectory,
#'     "hpc/fid/0_A20/1/1SRef/fid"), filterZeroIntensities=TRUE)
#'
#' ## plot spectrum
#' plot(noHpcSpec$spectrum$mass, noHpcSpec$spectrum$intensity, type="l",
#'      col="red", xlim=c(1296, 1300))
#' lines(hpcSpec$spectrum$mass, hpcSpec$spectrum$intensity, type="l",
#'       col="green", xlim=c(1296, 1300))
#' legend(x="topright", legend=c("no hpc", "hpc"), col=c("red", "green"), lwd=1)
#'
#' ## show difference between .hpc and original HPC
#' ## load mzXML generated by Bruker Daltonics CompassXport 1.3.5
#' ## you could do it like this:
#' #library("readMzXmlData")
#' #cpSpecHpcMzXml <- readMzXmlFile(file.path(exampleDirectory,
#' #  "hpc/mzXML/hpc.mzXML"))
#'
#' ## or easily use:
#' data(cpSpecHpcMzXml)
#'
#' ## reduce R double precision to single precision because our CompassXport 1.3.5
#' ## supports only mzXML with precision=32 (only for compatibility reasons)
#' noHpcSpec$spectrum$mass32 <-
#'  readBrukerFlexData:::.double2singlePrecision(noHpcSpec$spectrum$mass)
#' hpcSpec$spectrum$mass32 <-
#'  readBrukerFlexData:::.double2singlePrecision(hpcSpec$spectrum$mass)
#'
#' ## calculate deviance
#' d <- noHpcSpec$spectrum$mass32-cpSpecHpcMzXml$spectrum$mass
#' dHPC <- hpcSpec$spectrum$mass32-cpSpecHpcMzXml$spectrum$mass
#'
#' ## a little summary
#' cat("without .hpc:\n",
#'     "not matching: ", length(cpSpecHpcMzXml$spectrum$mass[d!=0]), " of ",
#'     length(cpSpecHpcMzXml$spectrum$mass), "; range: ",
#'     range(abs(d[d!=0])), "\nwith .hpc:\n",
#'     "not matching: ", length(cpSpecHpcMzXml$spectrum$mass[dHPC!=0]), " of ",
#'     length(cpSpecHpcMzXml$spectrum$mass), "; range: ",
#'     range(abs(d[dHPC!=0])), "\n")
#'
#' ##
#' ## doing things manually
#' ##
#' hpcMass <- readBrukerFlexData:::.hpc(mass=noHpcSpec$spectrum$mass,
#'  minMass=noHpcSpec$metaData$hpc$limits["minMass"],
#'  maxMass=noHpcSpec$metaData$hpc$limits["maxMass"],
#'  hpcCoefficients=noHpcSpec$metaData$hpc$coefficients)
#'
## TODO:
##  - internal calibration (or something like that)
##  - maybe need the following:
##      ##$Hpcgc0= 237.647814332814
##      ##$Hpcgc2= -0.0467016003168749
##
.hpc <- function(mass, minMass, maxMass, hpcCoefficients) {
  ## only defined for mass between minMass and maxMass
  ## QUESTION: Do we exclude values after or before calibration?
  ## ANSWER: I don't know how compassXport do it.
  ##         before: reduce calculation time
  ##         after: produce better results
  ##          I also don't know whether to use >= or > (<=, <).
  hpcRange <- mass >= minMass & mass <= maxMass
  m <- mass[hpcRange]

  ## correction=c[0] + c[1]*cal_mass^1 + ... + c[n]*cal_mass^n
  ## mass=cal_mass - correction
  l <- length(hpcCoefficients) - 1
  m <- sapply(m, function(x) {
    return(x - sum(hpcCoefficients * x^(0:l)))
  })

  mass[hpcRange] <- m

  return(mass)
}

