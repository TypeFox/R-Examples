#' Plateau test function for TL dating
#'
#' This function performs the plateau test for TL curves (Ln/Lx).
#'
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the TL curves used for the Plateau test.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'
#' @details
#' \bold{Plotting parameters} \cr
#' The plotting parameters are:  \cr
#' \describe{
#'  \item{\code{plot.Tmin}}{
#'    \link{numeric}: Lowest temperature plotted.}
#'  \item{\code{plot.Tmax}}{
#'    \link{numeric}: Highest temperature plotted.}
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#' See also \link{plot_TL.plateau}. \cr
#'
#' @return
#'  The results are plotted using \link{plot_TL.plateau}. \cr
#'
#' @references
#'  Aitken, M.J. (1985) Thermoluminescence Dating, Academic Press, London \cr
#'
#' @seealso
#'  \link{calc_TL.LxTx},
#'  \link{calc_TL.plateau},
#'  \link{analyse_TL.MAAD}
#'
#' @author David Strebler, University of Cologne (Germany)
#'
#' @export analyse_TL.plateau
#'
analyse_TL.plateau <- function(

  object,

  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA,
                           plateau.Tmin=0,
                           plateau.Tmax=NA,
                           no.plot=FALSE)
){

  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(object)){
    stop("[analyse_TL.SAR] Error: Input 'object' is missing.")

  }else if (!is(object,"TLum.Analysis")){
    stop("[analyse_TL.MAAD] Error: Input 'object' is not of type 'TLum.Analysis'.")
  }
  # ------------------------------------------------------------------------------

  sample.name <- as.character(object@records[[1]]@metadata$SAMPLE)

  data <- calc_TL.LxTx(object)

  temperatures  <- get_TLum.Results(data, "Temperatures")

  nPoints <- length(temperatures)
  Tmax <- max(temperatures)
  Tstep <- Tmax/nPoints

  dTypes <- get_TLum.Results(data, "Datatypes")
  doses <- get_TLum.Results(data, "Doses")
  Lx <- as.matrix(get_TLum.Results(data, "Lx"))
  Lx.error <- as.matrix(get_TLum.Results(data, "Lx.error"))
  LxTx <- as.matrix(get_TLum.Results(data, "LxTx"))
  LxTx.error <- as.matrix(get_TLum.Results(data, "LxTx.error"))

  # Additive curves
  temp.data <- calc_TL.MAAD.separate(doses=doses,
                                     Lx=Lx,
                                     Lx.error=Lx.error,
                                     dTypes=dTypes)

  aDoses <- get_TLum.Results(temp.data, "aDoses")
  aNames <- get_TLum.Results(temp.data, "aNames")
  aLx <- as.matrix(get_TLum.Results(temp.data, "aLx"))
  aLx.error <- as.matrix(get_TLum.Results(temp.data, "aLx.error"))

  temp.data <- calc_TL.MAAD.separate(doses=doses,
                                     Lx=LxTx,
                                     Lx.error=LxTx.error,
                                     dTypes=dTypes)

  aLxTx <- as.matrix(get_TLum.Results(temp.data, "aLx"))
  aLxTx.error <- as.matrix(get_TLum.Results(temp.data, "aLx.error"))

  # Averages
  temp.data <- calc_TL.MAAD.average(names=aNames,
                                    doses=aDoses,
                                    Lx=aLx,
                                    Lx.error=aLx.error)

  uDoses <- get_TLum.Results(temp.data, "doses")
  uNames <- get_TLum.Results(temp.data, "names")
  aLx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
  aLx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))

  temp.data <- calc_TL.MAAD.average(names=aNames,
                                    doses=aDoses,
                                    Lx=aLxTx,
                                    Lx.error=aLxTx.error)

  aLxTx.a <- as.matrix(get_TLum.Results(temp.data, "Lx"))
  aLxTx.a.error <- as.matrix(get_TLum.Results(temp.data, "Lx.error"))

  #Plateau test (Additive step)
  temp.data <- calc_TL.plateau(Ln=aLx.a[,uNames == "N"],
                               Ln.error=aLx.a.error[,uNames == "N"],
                               Lx=aLx.a[,uNames != "N" & uDoses != 0],
                               Lx.error=aLx.a.error[,uNames != "N"& uDoses != 0])

  aLx.plateau <- get_TLum.Results(temp.data,"LnLx")

  temp.data <- calc_TL.plateau(Ln=aLxTx.a[,uNames == "N"],
                               Ln.error=aLxTx.a.error[,uNames == "N"],
                               Lx=aLxTx.a[,uNames != "N" & uDoses != 0],
                               Lx.error=aLxTx.a.error[,uNames != "N"& uDoses != 0])

  aLxTx.plateau <- get_TLum.Results(temp.data,"LnLx")

  #----------------------------------------------------------------------------------------------------------------
  #Plot results
  #----------------------------------------------------------------------------------------------------------------

  no.plot <- plotting.parameters$no.plot

  # ------------------------------------------------------------------------------
  #Check values
  if(is.null(no.plot) || is.na(no.plot) || !is.logical(no.plot)){
    no.plot <- FALSE
  }
  # ------------------------------------------------------------------------------

  if(!no.plot){
    plot_TL.plateau(sample.name=sample.name,
                    temperatures=temperatures,
                    names=aNames,
                    doses=aDoses,
                    Lx=aLx,
                    Lx.a=aLx.a,
                    Lx.plateau=aLx.plateau,
                    LxTx=aLxTx,
                    LxTx.a=aLxTx.a,
                    LxTx.plateau=aLxTx.plateau,
                    plotting.parameters=plotting.parameters
    )
  }

  # ------------------------------------------------------------------------------

  new.TLum.Results.analyse_TL.plateau <- set_TLum.Results(data = list(names = uNames,
                                                                      doses = uDoses,
                                                                      Lx = aLx.a,
                                                                      Lx.error = aLx.a.error,
                                                                      Lx.plateau = aLx.plateau,
                                                                      LxTx=aLxTx,
                                                                      LxTx.a=aLxTx.a,
                                                                      LxTx.plateau=aLxTx.plateau)
  )

  return(new.TLum.Results.analyse_TL.plateau)
}




