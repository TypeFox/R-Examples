#' Aligning the TL peaks
#'
#' This function detects the peak position for each TL curve of the object and aligns them.
#' It uses the average of the testdose maximum positions as reference for the new peak position.
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the TL curves that have to be aligned.
#' @param aligning.parameters
#'  \link{list} (with default): list containing the aligning parameters. See details.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'
#'@details
#'  \bold{Aligning parameters} \cr
#'  The aligning parameters are:  \cr
#'  \describe{
#'  \item{\code{peak.Tmin}}{
#'    \link{numeric}: Lower boundary for looking at the peak maximum position.}
#'  \item{\code{peak.Tmax}}{
#'    \link{numeric}: Upper boundary for looking at the peak maximum position.}
#'  \item{\code{no.testdose}}{
#'    \link{logical}: If \code{TRUE}, the function will use the Lx curves rather than the Tx curves as reference for the peak maximum position.}
#'  }
#'
#'  \bold{Plotting parameters} \cr
#'  The plotting parameters are:  \cr
#'  \describe{
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#'
#' @return
#'  This function provides a new \code{\linkS4class{TLum.Analysis}} object with the same TL curves but aligned. \cr
#'  It also plots the original TL curves, the TL curves used to determine the peak maximum position, and the shiffted TL curves using \link{plot_align.peaks}.
#'
#' @seealso
#'  \link{plot_align.peaks}.
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export mod_align.peaks

mod_align.peaks <- function(

  object,

  aligning.parameters=list(peak.Tmin=0,
                           peak.Tmax=NA,
                           no.testdose=FALSE),

  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA,
                           no.plot=FALSE)

){

  C_TESTDOSE <- "Testdose"
  C_BACKGROUND <- "Background"

  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(object)){
    stop("[mod_align.peaks] Error: Input 'object' is missing.")
  }else if (!is(object,"TLum.Analysis")){
    stop("[mod_align.peaks] Error: Input 'object' is not of type 'TLum.Analysis'.")
  }

  if(!is.list(aligning.parameters)){
    stop("[mod_align.peaks] Error: Input 'aligning.parameters' is not of type 'list'.")
  }

  if(!is.list(plotting.parameters)){
    stop("[mod_align.peaks] Error: Input 'plotting.parameters' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  peak.Tmin <- aligning.parameters$peak.Tmin
  peak.Tmax <- aligning.parameters$peak.Tmax
  no.testdose <- aligning.parameters$no.testdose

  no.plot <- plotting.parameters$no.plot


  protocol <- object@protocol

  nRecords <- length(object@records)

  nChannel <- vector()

  for(i in 1: nRecords){
    temp.record <- object@records[[i]]

    temp.nChannel <- length(temp.record@data)
    nChannel <- c(nChannel, temp.nChannel)
  }

  nChannel <- unique(nChannel)


  temperatures <- object@records[[1]]@temperatures

  Tmax <- max(temperatures)

  # ------------------------------------------------------------------------------
  # Value check
  if(length(nChannel) != 1){
    stop("[mod_align.peaks] Error: All records do not have the same number of channels.")
  }

  for (i in 1: nRecords){
    temp.record <- object@records[[i]]

    temp.temperatures <- temp.record@temperatures

    if(!identical(temperatures, temp.temperatures)){
      stop("[mod_align.peaks] Error: All records do not have the same Temperatures vector.")
    }
  }

  if(!is.numeric(peak.Tmin) || !is.finite(peak.Tmin)){
    peak.Tmin <- 0
  }else if(peak.Tmin < 0){
    stop("[mod_align.peaks] Error: peak.Tmin < 0.")
  }else if(peak.Tmin > Tmax){
    stop("[mod_align.peaks] Error: peak.Tmin > Tmax.")
  }else if(peak.Tmin > peak.Tmax){
    stop("[mod_align.peaks] Error: peak.Tmin > peak.Tmax")
  }

  if(!is.numeric(peak.Tmax) || !is.finite(peak.Tmax)){
    peak.Tmax <- Tmax
  }else if(peak.Tmax > Tmax){
    peak.Tmax <- Tmax
  }

  if(is.null(no.testdose) || is.na(no.testdose) || !is.logical(no.testdose)){
    no.testdose <- FALSE
  }

  if(is.null(no.plot) || is.na(no.plot) || !is.logical(no.plot)){
    no.plot <- FALSE
  }
  # ------------------------------------------------------------------------------

  min <- 1
  max <- length(temperatures)

  for(i in 1:length(temperatures)){
    if(temperatures[i]<=peak.Tmin){
      min <- i
    }
    if(temperatures[i]<=peak.Tmax){
      max <- i
      }
  }

  #----------------------------------------------------------------------------------------------
  #Separation Tx and Lx
  #----------------------------------------------------------------------------------------------

  Lx <- vector()
  Lx.error <- vector()
  Tx <- vector()
  Tx.error <- vector()
  TL <- vector()
  TL.error <- vector()

  for(i in 1:nRecords){
    temp.record <- object@records[[i]]

    temp.dtype <- temp.record@metadata$DTYPE

    temp.curve <- temp.record@data
    temp.curve.error <- temp.record@error

    TL <- cbind(TL,temp.curve)
    TL.error <- cbind(TL.error,temp.curve.error)

    if(temp.dtype != C_TESTDOSE) {
      Lx <- cbind(Lx,temp.curve)
      Lx.error <- cbind(Lx.error,temp.curve.error)
    }else{
      Tx <- cbind(Tx,temp.curve)
      Tx.error <- cbind(Tx.error,temp.curve.error)
    }
  }

  # ---------------------------------------------------
  # Values check
  if(length(Tx) == 0 && !no.testdose){
    warning("[mod_align.peaks] Warning: no testdoses registered.")
  }else if(length(Tx) != length(Lx)){
    warning("[mod_align.peaks] Warning: The Lx and Tx matrix have different size.")
  }

  # --------------------------------------------------
  #Reference peak

  if(length(Tx) > 0 || !no.testdose){
    Tx.a <- rowMeans(Tx)
    Tx.a.smooth <- smooth.spline(Tx.a)
    Tx.a.lim <- Tx.a.smooth$y[min:max]
    max.Tx.lim <- which.max(Tx.a.lim)

    new.peak <-  max.Tx.lim + (min-1)

  }else{
    Lx.a <- rowMeans(Lx)
    Lx.a.smooth <- smooth.spline(Lx.a)
    Lx.a.lim <- Lx.a.smooth$y[min:max]
    max.Lx.lim <- which.max(Lx.a.lim)

    new.peak <-  max.Lx.lim + (min-1)
  }

  Tpeak <- temperatures[new.peak]


  #----------------------------------------------------------------------------------------------
  #Peak Shift
  #----------------------------------------------------------------------------------------------

  new.TL <- vector()
  new.TL.error <- vector()

  for(i in 1:nRecords){

    temp.record <- object@records[[i]]
    temp.dtype <- temp.record@metadata$DTYPE

    temp.curve <- TL[,i]
    temp.curve.error <- TL.error[,i]

    #Shift Calculation
    temp.smooth <- smooth.spline(temp.curve)
    temp.smooth.lim <- temp.smooth$y[min:max]

    max.temp.smooth.lim <- which.max(temp.smooth.lim)
    temp.max <- max.temp.smooth.lim + (min-1)

    if(temp.curve[temp.max]> 120 && temp.dtype != C_BACKGROUND){
      temp.shift <- temp.max - new.peak
    }else{
      temp.shift <- 0     # Do not shift BG signal or a signal with very low intensity (<120 counts/channel)
    }


    #Peak shifting
    temp.new.curve <- integer()
    temp.new.curve.error <- integer()

    for(j in 1:length(temp.curve)){
      k <- j + temp.shift

      if(k > 0 && k <= length(temp.curve)){
        temp.new.curve[j] <- temp.curve[k]
        temp.new.curve.error[j] <- temp.curve.error[k]
      }else{
        temp.new.curve[j] <- NA
        temp.new.curve.error[j] <- NA
      }
    }

    new.TL <- cbind(new.TL,temp.new.curve)
    new.TL.error <- cbind(new.TL.error,temp.new.curve.error)
  }

  #----------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------

  if(!no.plot){
    if(length(Tx) > 0 || !no.testdose){
      plot_align.peaks(temperatures=temperatures,
                       old.TL=TL,
                       new.TL=new.TL,
                       ref.TL=Tx,
                       pos.peak=Tpeak,
                       plotting.parameters=plotting.parameters)
    }else{
      plot_align.peaks(temperatures=temperatures,
                       old.TL=TL,
                       new.TL=new.TL,
                       ref.TL=Lx,
                       pos.peak=Tpeak,
                       plotting.parameters=plotting.parameters)
    }
  }

  #----------------------------------------------------------------------------------------------
  # New TLum.Analysis
  #----------------------------------------------------------------------------------------------

  new.records <- list()

  for(i in 1:nRecords){

    temp.record <- object@records[[i]]

    temp.record@data <- new.TL[,i]
    temp.record@error <- new.TL.error[,i]

    new.records <- c(new.records, temp.record)

  }

  new.analysis <- set_TLum.Analysis(records = new.records,
                                    protocol = protocol)

  return(new.analysis)
}
