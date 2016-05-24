#' Extract TL
#'
#' This function provides a new \code{\linkS4class{TLum.Analysis}} object containing only the TL curves.
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the initial TL curves.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#' @param record.parameters
#'  \link{list} (with default): list containing the record parameters. See details.
#'
#'@details
#'  \bold{Plotting parameters} \cr
#'  The plotting parameters are:  \cr
#'  \describe{
#'    \item{\code{no.plot}}{
#'      \link{logical}: If \code{TRUE}, the results will not be plotted.}
#'  }
#'
#'  \bold{Record parameters} \cr
#'  The record parameters are:  \cr
#'  \describe{
#'    \item{\code{includePreheat}}{
#'      \link{logical}: If \code{TRUE}, the preheat was included in the TL recording. If \code{FALSE}, the preheat was recorded separately.}
#'    \item{\code{recDuringPreheatRamp}}{
#'      \link{logical}: Only used when \code{includePreheat} is \code{TRUE}. If \code{TRUE}, the signal was recorded during the preheat ramp.}
#'      \item{\code{recDuringPreheatPlateau}}{
#'      \link{logical}: Only used when \code{includePreheat} is \code{TRUE}. If \code{TRUE}, the signal was recorded during the preheat plateau.}
#'  }
#'
#' @return
#'  This function provides a new \code{\linkS4class{TLum.Analysis}} with only the TL curve. \cr
#'  It also plots the TL curves using \link{plot_extract.TL}.
#'
#' @seealso
#'  \link{plot_extract.TL}
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export mod_extract.TL

mod_extract.TL <- function(

  object,

  plotting.parameters=list(no.plot=FALSE),

  record.parameters=list(separatePreheat=TRUE,
                         recDuringPreheatRamp=FALSE,
                         recDuringPreheatPlateau=FALSE)

){
  C_TL <- "TL"

  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(object)){
    stop("[mod_extract.TL] Error: Input 'object' is missing.")
  }else if (!is(object,"TLum.Analysis")){
    stop("[mod_align.peaks] Error: Input 'object' is not of type 'TLum.Analysis'.")
  }

  if(!is.list(plotting.parameters)){
    stop("[mod_extract.TL] Error: Input 'plotting.parameters' is not of type 'list'.")
  }

  if(!is.list(record.parameters)){
    stop("[mod_extract.TL] Error: Input 'record.parameters' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  no.plot <- plotting.parameters$no.plot

  # ------------------------------------------------------------------------------
  # Value check
  if(is.null(no.plot) || is.na(no.plot) || !is.logical(no.plot)){
    no.plot <- FALSE
  }
  # ------------------------------------------------------------------------------

  new.protocol <- object@protocol
  records <- object@records

  nRecords <- length(records)

  new.records <- list()
  rejected.records <- list()

  TL <- list()
  temperatures <- list()

  test.TL <- vector()

  new.id <- 1

  for(i in 1:nRecords){
    temp.curve <- records[[i]]

    temp.metadata <- temp.curve@metadata

    temp.ltype <- temp.metadata$LTYPE

    temp.data <- list(temp.curve@data)
    temp.temperatures <- list(temp.curve@temperatures)

    if(temp.ltype == C_TL){
      temp.test <- TRUE

      TL <- c(TL, temp.data)
      temperatures <- c(temperatures, temp.temperatures)

      new.curve <- temp.curve
      new.curve@metadata$ID <- new.id
      new.id <- new.id+1

      new.records <- c(new.records, new.curve)
    }else{
      temp.test <- FALSE
      rejected.records <- c(rejected.records, temp.curve)
    }

    test.TL <- c(test.TL,temp.test)
  }

  new.TLum.Analysis <- set_TLum.Analysis(records= new.records,
                                         protocol=new.protocol)

  #--------------------------------------------------------------------------------------------------------
  #Plot results
  #--------------------------------------------------------------------------------------------------------

  if(!no.plot){
    plot_extract.TL(TL=TL, temperatures=temperatures)
  }

  # -------------------------------------------------------------------------------------------------------------------
  #New TLum.Analyis
  # -------------------------------------------------------------------------------------------------------------------


  return(new.TLum.Analysis)

}
