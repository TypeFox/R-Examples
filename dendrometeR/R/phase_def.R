#' Define stem-cyclic phases
#'
#' @description The function identifies and assigns each timestamp to one of the three distinct phases of contraction, expansion and stem-radius increment (Deslauriers et al. 2011) for dendrometer series from a \code{data.frame} with gap-free dendrometer data.
#'
#' @usage phase_def(dm.gpf, resolution = dendro.resolution(dm.gpf),
#'           shapeSensitivity = 0.6, minmaxDist = 0.2, minmaxSD = 2,
#'           radialIncrease = "max")
#'
#' @param dm.gpf a \code{data.frame} with either gap-free or gap-filled dendrometer series as produced by \code{\link{fill_gaps}}.
#' @param resolution a \code{numeric} specifying the resolution of the dendrometer data in seconds. Defaults to the resolution of \code{dm.gpf} as calculated using \code{\link{dendro.resolution}}.
#' @param shapeSensitivity a \code{numeric} specifying a time window, defined as proportion of a single day. Within this time window possible extrema points (i.e. minimum and maximum) in dendrometer measurements are searched for. Defaults to 0.6 (60\% of a day). See details for further explanation.
#' @param minmaxDist a \code{numeric} specifying the minimum temporal distance between consecutive minimum and maximum points (i.e. in the x direction). Defaults to 0.2 (20\% of a day).
#' @param minmaxSD a \code{numeric} specifying the minimum difference between consecutive minimum and maximum points expressed as a number of standard deviations (i.e. in the y direction). Defaults to 2.
#' @param radialIncrease a \code{character} string of \code{"max", "min", "mid"}, specifying when the stem-radius increment phase should start, with \code{"max"} as the most, and \code{"min"} as the least conservative approach; \code{"mid"} is in between. See details for further explanation.
#'
#' @details The function defines the stem-cyclic phases of contraction, expansion, and stem-radius increment, as described in Deslauriers et al. (2011). The function is a more robust version of the original SAS routine, as its architecture allows to handle noisy data as well.
#'
#' First, the function searches for minimum and maximum points within a daily time window as specified by \code{shapeSensitivity}. Then, the original dendrometer series are offset by \code{(1 - shapeSensitivity) / 2} in both directions to assure whether the identified extrema are indeed the extrema of cyclic phases. A comparison between the original and offset series allows to select all appropriate minimum and maximum values.
#'
#' The arguments \code{minmaxDist} and \code{minmaxSD} specify the temporal distance and the minimum difference between consecutive minimum and maximum points (i.e. in x and y direction), respectively. The argument \code{radialIncrease} determines from which moment on data points should be assigned to the stem-radius increment phase: when points are continuously above the previous maximum (\code{"max"}), when a single data point is above the previous maximum (\code{"min"}), or right in between \code{"min"} and \code{"max"} (\code{"mid"}).
#'
#' @return
#' The function returns a \code{data.frame} with numbers indicating the different stem-cyclic phases: (1) contraction, (2) expansion, (3) stem-radius increment for each timestamp.
#'
#' @author Marko Smiljanic
#'
#' @references Deslauriers, A., Rossi, S., Turcotte, A., Morin, H. and Krause, C. (2011) A three-step procedure in SAS to analyze the time series from automatic dendrometers. \emph{Dendrochronologia} 29: 151-161.
#'
#' @examples
#' data(dmCD)
#' dm.phase <- phase_def(dmCD)
#'
#' @export phase_def
#'
phase_def <- function(dm.gpf, resolution = dendro.resolution(dm.gpf), shapeSensitivity = 0.6, minmaxDist = 0.2, minmaxSD = 2, radialIncrease = "max")
{
  ### Recursive call if there is more than one dendrometer series ###
  if(class(dm.gpf) != "numeric") {
    dm.gpf <- as.data.frame(dm.gpf)
    output <- as.data.frame(matrix(nrow = nrow(dm.gpf), ncol = ncol(dm.gpf)))

    if(is.dendro(dm.gpf)) {
      for(d in 1:ncol(dm.gpf)) {
        den <- dm.gpf[[d]]
        res <- phase_def(den,resolution,shapeSensitivity, minmaxDist, minmaxSD, radialIncrease)
        rowRange <- (nrow(output)-length(res)+1):nrow(output)
        output[rowRange, d] <- res
      }
      names(output) <- names(dm.gpf)
      rownames(output) <- rownames(dm.gpf)
    }
    else {
      stop("there is a problem with the input data")
    }
    return(output)
  }
  if(any(is.na(dm.gpf))) {
    stop("NA values detected in dendrometer series. Run fill_gaps first")
  }
  ### Detect local minima and maxima ###
  localExtrema <- function(dendro = parent.frame()$dm.gpf) {
    minima <- vector()
    maxima <- vector()
    i <- 2
    while(i < length(dendro)) {
      dat <- dendro[i]
      datAfter <- dendro[i+1]
      datBefore <- dendro[i-1]
      if(dat == datAfter) {
        dat <- dendro[i] <- dendro[i] + (dat-datBefore)/2
      }
      if(dat == datBefore) {
        iStart <- i
        while(dendro[i] == dat) {
          i <- i+1
        }
        pos <- ((iStart-1):i)
        vals <- (pos-1)*(dendro[i]-dat)/(i-iStart+1)+dat
      }
      if(!is.na(datAfter) && !is.na(datBefore)) {
        if(dat > datAfter && dat > datBefore) {
          maxima <- c(maxima, i)
        }
        else if(dat < datAfter && dat < datBefore) {
          minima <- c(minima, i)
        }
      }
      else {
        stop("NA values detected in dendrometer series")
      }
      i <- i+1
    }

    result <- list()
    result$minima <- minima
    result$maxima <- maxima
    return(result)

  }

  ### select minima/maxima points for every day depending on the resolution ###
  extrema <- function(dendro = parent.frame()$dm.gpf, minima = parent.frame()$minima, maxima = parent.frame()$maxima, shapeSensitivity = parent.frame()$shapeSensitivity) {
    dayLen <- 86400/resolution
    dayNumb <- round(length(dendro)/dayLen)
    dayShape <- vector()
    minimaNew <- vector()
    maximaNew <- vector()
    for(day in 1:dayNumb) {
      dayStart <- (day-1)*dayLen+1
      dayEnd <- min(dayStart+dayLen, length(dendro))
      dendroDay <- dendro[dayStart:dayEnd]

      minDendroDay <- min(dendroDay, na.rm=TRUE)
      minDendroDay <- which(dendroDay == minDendroDay)

      maxDendroDay <- max(dendroDay, na.rm=TRUE)
      maxDendroDay <- which(dendroDay == maxDendroDay)

      lowSens <- dayLen*shapeSensitivity
      highSens <- (1-shapeSensitivity)*dayLen
      if(minDendroDay > lowSens && minDendroDay < highSens) {
        dayShp <- 0
      }
      else if(maxDendroDay > lowSens && maxDendroDay < highSens) {
        dayShp <- 1
      }
      else {
        dayShp <- NA
      }
      dayShape <- c(dayShape, dayShp)

      if(!is.na(dayShp)) {
        dendroMinimaID <- minDendroDay+dayStart-1
        dendroMaximaID <- maxDendroDay+dayStart-1
        dendroMinimaID <- dendroMinimaID[[1]]
        dendroMaximaID <- dendroMaximaID[[1]]
        minimaNew <- c(minimaNew, minima[which(minima==dendroMinimaID)])
        maximaNew <- c(maximaNew, maxima[which(maxima==dendroMaximaID)])
      }
    }
    result <- list(minimaNew,maximaNew,dayShape)
    names(result) <- c("minimaNew", "maximaNew", "dayShape")
    return(result)
  }

  shapeSensitivity = (1-shapeSensitivity)/2
  locExt <- localExtrema()
  minima <- locExt$minima
  maxima <- locExt$maxima
  extVals <- extrema()

  # offset original series depending on the shapeSensitivity variable and redo the extrema detection
  # TODO: Offset might be too big currently. Needs further testing
  dayLen <- 86400/resolution
  locExtOffset <- localExtrema(dm.gpf[floor(dayLen*shapeSensitivity):length(dm.gpf)])
  extValsOffset <- extrema(dm.gpf[floor(dayLen*shapeSensitivity):length(dm.gpf)], locExtOffset$minima, locExtOffset$maxima)

  ### Depending on the number of indeterminate day curve shapes select either original or offset series
  origShape <- extVals$dayShape
  offsetShape <- extValsOffset$dayShape
  origShapeNAs <- sum(is.na(origShape))
  offsetShapeNAs <- sum(is.na(offsetShape))

  lenConc <- length(origShape) - length(which(origShape == 0))
  lenConv <- length(origShape) - length(which(origShape == 1))
  lenConcOff <- length(offsetShape) - length(which(offsetShape == 0))
  lenConvOff <- length(offsetShape) - length(which(offsetShape == 1))
  diffCC <- abs(lenConc - lenConv)
  diffCCOff <- abs(lenConcOff - lenConvOff)

  dm.orig <- dm.gpf
  if(offsetShapeNAs < origShapeNAs) {
    dm.gpf <- dm.gpf[floor(dayLen*shapeSensitivity):length(dm.gpf)]
    minima <- locExtOffset$minima
    maxima <- locExtOffset$maxima
    extVals <- extValsOffset
  }
  else if(offsetShapeNAs == origShapeNAs) {
    if(diffCCOff > diffCC) {
      dm.gpf <- dm.gpf[floor(dayLen*shapeSensitivity):length(dm.gpf)]
      minima <- locExtOffset$minima
      maxima <- locExtOffset$maxima
      extVals <- extValsOffset
    }
  }

  if(length(dm.orig) > length(dm.gpf)) {
    reps <- length(dm.orig)-length(dm.gpf)
    dm.gpf <- c(rep(NA, reps), dm.gpf)
    extVals$minimaNew <- extVals$minimaNew+reps
    extVals$maximaNew <- extVals$maximaNew+reps
    minima <- minima+reps
    maxima <- maxima+reps
  }

  minimaNew <- extVals$minimaNew
  maximaNew <- extVals$maximaNew
  dayShape <- extVals$dayShape

  #	par(ask=FALSE) #DEBUG
  #	plot(dm.gpf, type="l") #DEBUG

  ### Depending on the ratio of convex to concave day curve shapes maxima or minima will hold true values.
  ### Fill the other vector from its local minima values
  lenConc <- length(dayShape) - length(which(dayShape == 0))
  lenConv <- length(dayShape) - length(which(dayShape == 1))
  ratioCC <- lenConc/lenConv

  if(ratioCC <= 1) {
    maximaNew <- vector()
    for(mini in 1:(length(minimaNew)-1)) {
      minimaStart <- minimaNew[mini]
      minimaEnd <- minimaNew[mini+1]
      maximas <- maxima[which(maxima > minimaStart & maxima < minimaEnd)]
      dendroMax <- dm.gpf[maximas]
      maximaNew <- c(maximaNew, maximas[which(dendroMax == max(dendroMax))][1])
    }
  }
  else if(ratioCC > 1) {
    minimaNew <- vector()
    for(maxi in 1:(length(maximaNew)-1)) {
      maximaStart <- maximaNew[maxi]
      maximaEnd <- maximaNew[maxi+1]
      minimas <- minima[which(minima > maximaStart & minima < maximaEnd)]
      dendroMin <- dm.gpf[minimas]
      minimaNew <- c(minimaNew, minimas[which(dendroMin == min(dendroMin))][1])
    }
  }

  ### if the minima maxima pair is too close to each other delete both
  for(i in 1:length(minimaNew)) {
    mini <- minimaNew[i]
    potMaxiID <- which(abs(maximaNew-mini) < dayLen*minmaxDist)
    if(length(potMaxiID) != 0 ) {
      potMaxi <- maximaNew[potMaxiID]
      if(mini != length(dm.gpf) && potMaxi != length(dm.gpf) && mini != 1 && potMaxi != 1) {
        minimaNew[i] <- maximaNew[potMaxiID] <- NA
        minimaNew <- na.omit(minimaNew)
        maximaNew <- na.omit(maximaNew)
      }
    }
  }

  ### if the difference between the minima and maxima in minima maxima pair is too small delete both
  extremaPoints <- as.data.frame(rbind(cbind(maximaNew,0),cbind(minimaNew, 1)))
  names(extremaPoints) <- c("Position", "type")
  extremaPoints <- extremaPoints[order(extremaPoints$Position),]
  extremaValues <- as.data.frame(cbind(dm.gpf[extremaPoints$Position], extremaPoints$type))
  extremaDiff <- abs(diff(extremaValues[,1]))
  extremaSD <- (extremaDiff-mean(extremaDiff))/sd(extremaDiff)
  extremaPot <- which(extremaSD < (-minmaxSD))
  for(ext in extremaPot) {
    extPot <- extremaPoints[ext:(ext+1),]
    if(diff(extPot$type) == 1) {
      for(p in extPot$Position){
        maximaNew[which(maximaNew == p)] <- NA
        minimaNew[which(minimaNew == p)] <- NA
      }
    }
  }
  minimaNew <- na.omit(minimaNew)
  maximaNew <- na.omit(maximaNew)
  #points(y=dm.gpf[minimaNew],x=minimaNew, col=2) #DEBUG
  #points(x=maximaNew, y=dm.gpf[maximaNew], col=3) #DEBUG

  ### Assign the growth phases ###
  minima <- minimaNew
  maxima <- maximaNew
  minimaLen <- length(minima)
  maximaLen <- length(maxima)
  minimaFirst <- minima[1] < maxima[1]
  maximaLast <- minima[length(minima)] < maxima[length(maxima)]
  phase <- vector()

  # Depending on the order which comes first in dendrometer series (minima or maxima)
  # this if statement goes through extrema points and assigns the phases to vectors
  # currently in the if statement one block of code is copied and minima replaced with maxima
  # TODO: Avoid this copy with minima or maxima being a reference. If that does not work use a function
  if(minimaFirst) {
    for(i in 1:length(minima)) {
      if(!is.na(maxima[i])) {
        phase[(minima[i]+1):maxima[i]] <- 2
        if(i > 1) {
          ref <- dm.gpf[maxima[i-1]]
          dotsPhase3 <- which(dm.gpf[minima[i]:maxima[i]] > ref)
          phaseStart <- NULL
          if(length(dotsPhase3) > 0 && radialIncrease == "min") {
            phaseStart <- min(dotsPhase3)
          }
          else if(length(dotsPhase3) > 0 && radialIncrease == "max") {
            if(length(dotsPhase3) ==  1) {
              phaseStart <- dotsPhase3
            }
            else {
              while(length(unique(diff(dotsPhase3))) != 1) {
                dotsPhase3 <- dotsPhase3[2:length(dotsPhase3)]
              }
              if(unique(diff(dotsPhase3)) == 1) {
                phaseStart <- min(dotsPhase3)
              }
              else {
                phaseStart <- max(dotsPhase3)
              }
            }
          }
          else if(length(dotsPhase3) > 0 && radialIncrease == "mid") {
            phaseMin <- min(dotsPhase3)
            if(length(dotsPhase3) ==  1) {
              phaseMax <- dotsPhase3
            }
            else {
              while(length(unique(diff(dotsPhase3))) != 1) {
                dotsPhase3 <- dotsPhase3[2:length(dotsPhase3)]
              }
              if(unique(diff(dotsPhase3)) == 1) {
                phaseMax <- min(dotsPhase3)
              }
              else {
                phaseMax <- max(dotsPhase3)
              }
            }
            phaseStart <- round(mean(c(phaseMin, phaseMax)))
          }
          if(!is.null(phaseStart)) {
            phaseStart <- phaseStart + minima[i] - 1
            phase[phaseStart:maxima[i]] <- 3
          }
        }
      }
      if(!is.na(minima[i+1]) && !is.na(maxima[i])) {
        phase[(maxima[i]+1):minima[i+1]] <- 1
      }
    }
  }
  else {
    if(maxima[1] < 1) {
      phase[1:maxima[1]] <- NA
    }
    for(i in 1:length(maxima)) {
      if(!is.na(minima[i])) {
        phase[(maxima[i]+1):minima[i]] <- 1
      }
      if(!is.na(maxima[i+1])) {
        phase[(minima[i]+1):maxima[i+1]] <- 2
        ref <- dm.gpf[maxima[i]]
        dotsPhase3 <- which(dm.gpf[minima[i]:maxima[i+1]] > ref)
        phaseStart <- NULL
        if(length(dotsPhase3) > 0 && radialIncrease == "min") {
          phaseStart <- min(dotsPhase3)
        }
        else if(length(dotsPhase3) > 0 && radialIncrease == "max") {
          if(length(dotsPhase3) == 1) {
            phaseStart <- dotsPhase3
          }
          else {
            while(length(unique(diff(dotsPhase3))) != 1) {
              dotsPhase3 <- dotsPhase3[2:length(dotsPhase3)]
            }
            if(unique(diff(dotsPhase3)) == 1) {
              phaseStart <- min(dotsPhase3)
            }
            else {
              phaseStart <- max(dotsPhase3)
            }
          }
        }
        else if(length(dotsPhase3) > 0 && radialIncrease == "mid") {
          phaseMin <- min(dotsPhase3)
          if(length(dotsPhase3) ==  1) {
            phaseMax <- dotsPhase3
          }
          else {
            while(length(unique(diff(dotsPhase3))) != 1) {
              dotsPhase3 <- dotsPhase3[2:length(dotsPhase3)]
            }
            if(unique(diff(dotsPhase3)) == 1) {
              phaseMax <- min(dotsPhase3)
            }
            else {
              phaseMax <- max(dotsPhase3)
            }
          }
          phaseStart <- round(mean(c(phaseMin, phaseMax)))
        }
        if(!is.null(phaseStart)) {
          phaseStart <- phaseStart + minima[i] - 1
          phase[phaseStart:maxima[i+1]] <- 3
        }
      }
    }
  }

  if(length(phase) != length(dm.gpf)) {
    phase[(length(phase)+1):length(dm.gpf)] <- NA
  }
  #DEBUG: points(x=1:length(phase), y=dendro[1:length(phase)], pch=19, col=phase, cex=1)
  return(phase)
}
