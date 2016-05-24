#' Creates a new TLum.Analysis object where the background was removed from the signal.
#'
#'
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the initial TL curves.
#' @param keep.background
#'  \link{logical} (with default): Parameter indicating if the background curve have to be kept or suppressed.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'
#'@details
#'  \bold{Plotting parameters} \cr
#'  The plotting parameters are:  \cr
#'  \describe{
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#'
#' @seealso
#'  \link{plot_substract.background}
#'
#' @return
#'  This function provides a new \code{\linkS4class{TLum.Analysis}} object with the TL curves after background subtraction. \cr
#'  It also plots the TL curves, the background curves and the background substracted curves using \link{plot_remove.preheat}.
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export mod_substract.background

mod_substract.background <- function(

  object,

  keep.background=FALSE,

  plotting.parameters=list(no.plot=FALSE)
){

  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(object)){
    stop("[mod_substract.background] Error: Input 'object' is missing.")
  }else if (!is(object,"TLum.Analysis")){
    stop("[mod_substract.background] Error: Input 'object' is not of type 'TLum.Analysis'.")
  }

  if(!is.logical(keep.background) || is.na(keep.background)){
    stop("[mod_substract.background] Error: Input 'keep.background' is not of type 'logical'.")
  }

  if(!is.list(plotting.parameters)){
    stop("[mod_substract.background] Error: Input 'plotting.parameters' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  protocol <- object@protocol

  nRecords <- length(object@records)

  no.plot <- plotting.parameters$no.plot

  # ------------------------------------------------------------------------------
  # Value check
  if(is.null(no.plot) || is.na(no.plot) || !is.logical(no.plot)){
    no.plot <- FALSE
  }
  # ------------------------------------------------------------------------------

  #Extract BG & TL
  test.background <- logical()

  BG <- vector()
  BG.error <- vector()
  BG.temperature <- vector()

  TL <- vector()
  TL.error <- vector()
  TL.temperature <- vector()


  for(i in 1:nRecords){
    temp.record <- object@records[[i]]

    temp.curve <- temp.record@data
    temp.curve.error <- temp.record@error

    temp.temperatures <- temp.record@temperatures

    temp.metadata <- temp.record@metadata
    temp.dtype <- temp.metadata$DTYPE

    if(temp.dtype == "Background"){

      test.background[i] <- TRUE

      BG <- cbind(BG,temp.curve)
      BG.error <- cbind(BG.error,temp.curve.error)
      BG.temperature <- cbind(BG.temperature, temp.temperatures)

    }else{
      test.background[i] <- FALSE

      TL <- cbind(TL,temp.curve)
      TL.error <- cbind(TL.error,temp.curve.error)
      TL.temperature <- cbind(TL.temperature, temp.temperatures)
    }
  }

  #----------------------------------------------------------------------------------------------
  #Background substraction
  #----------------------------------------------------------------------------------------------

  if(identical(TL.temperature,BG.temperature)){
    new.TL <- TL - BG
    new.TL.error <- sqrt(TL.error^2 + BG.error^2)
  }else{
    stop("[mod_substract.background] Error: TL & BG Temperature matrix do not match.")
  }

  temperatures <- TL.temperature[,1]
  for(i in 1:ncol(TL.temperature)){
    if(!identical(temperatures,TL.temperature[,i])){
      stop("[mod_substract.background] Error: All TL do not have the same temperature vector.")
    }
  }
  #----------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------

  if(no.plot==FALSE){
    plot_substract.background(old.TL=TL,
                              BG=BG,
                              new.TL=new.TL,
                              temperatures=temperatures)
  }


  #----------------------------------------------------------------------------------------------
  # New BinFileData
  #----------------------------------------------------------------------------------------------


  new.records <- list()

  if(keep.background == FALSE){

    temp.id <- 0

    for(i in 1:nRecords){
      temp.record <- object@records[[i]]

      if(test.background[i] == FALSE) {

        temp.id <- temp.id+1
        temp.record@metadata$ID <- temp.id

        temp.record@data <- new.TL[,temp.id]
        temp.record@error <- new.TL.error[,temp.id]

        new.records <- c(new.records,temp.record)
      }
    }
  }else{
    #If keep.background == TRUE... for "analyse_SAR.TL"
    for(i in 1:nRecords){
      temp.record <- object@records[[i]]

      if(test.background[i] == FALSE) {

        temp.record@data <- new.TL[,temp.id]
        temp.record@error <- new.TL.error[,temp.id]
      }

      new.records <- c(new.records,temp.record)
    }
  }

  # new Analysis
  new.TLum.Analysis <- set_TLum.Analysis(records= new.records,
                                         protocol=protocol)

  return(new.TLum.Analysis)
}
