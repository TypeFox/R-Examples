#' Script for the MAAD protocol
#'
#' This script calls a series of data functions before estimating the ED using the MAAD protocol.
#' It just requires the name of the files with the TL curves, the relative error on the measurement and the temperature boundaries for the signal integration.
#'
#' @param file.name
#'  \link{character} (\bold{required}): Name of the file containing the luminescence data.
#' @param eval.Tmin
#'  \link{integer} (\bold{required}): Temperature ( °C) of the lower boundary for the signal integration.
#' @param eval.Tmax
#'  \link{integer} (\bold{required}): Temperature (°C) of the upper boundary for the signal integration.
#' @param file.parameters
#'  \link{list} (with default): list containing the input/output parameters. See details.
#' @param aligning.parameters
#'  \link{list} (with default): list containing the aligning parameters. See details.
#' @param fitting.parameters
#'  \link{list} (with default): list containing the fitting parameters. See details.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#' @param rejection.criteria
#'  \link{list} (with default): list containing the rejection criteria (in \%). See details.
#' @param relative.error
#'  \link{numeric} (with default): Relative error of the TL signals.
#' @param remove.discs
#'  \link{numeric}  (with default): list containing the position of the aliquots to shall be removed.
#'
#' @details
#' \bold{File parameters} \cr
#' The file parameters are:  \cr
#' \describe{
#'  \item{\code{file.extension}}{
#'    \link{character} (with default): extension of the file containing the luminescence data (.bin or .binx)}
#'  \item{\code{folder.in}}{
#'    \link{character} (with default): Folder containing the file with the luminescene data.}
#'  \item{\code{folder.out}}{
#'    \link{character} (with default): Folder containing the file with the new luminescene data.}
#' }
#' see also \link{script_TL.pretreatment}.
#'
#'  \bold{Aligning parameters} \cr
#'  The aligning parameters are:  \cr
#'  \describe{
#'  \item{\code{peak.Tmin}}{
#'    \link{numeric}: Lower boundary for looking for the peak maximum position.}
#'  \item{\code{peak.Tmax}}{
#'    \link{numeric}: Upper boundary for looking for the peak maximum position.}
#'  \item{\code{no.testdose}}{
#'    \link{logical}: If \code{TRUE}, the function will use the Lx curves rather than the Tx curves as reference for the peak maximum position.}
#'  }
#'
#' \bold{Fitting parameters} \cr
#' The fitting parameters are:  \cr
#' \describe{
#'  \item{\code{method}}{
#'    \link{character}: Fitting method (\code{LIN}, \code{EXP}, \code{EXP+LIN} or \code{EXP+EXP}).}
#'  \item{\code{fit.weighted}}{
#'    \link{logical}: If the fitting is weighted or not.}
#'  \item{\code{fit.use.slope}}{
#'    \link{logical}: If the slope of the Q growth curve is reused for the sublinearity correction.}
#'  \item{\code{fit.aDoses.min}}{
#'    \link{numeric}: Lowest additive dose used for the fitting.}
#'  \item{\code{fit.aDoses.max}}{
#'    \link{numeric}: Highest additive dose used for the fitting.}
#'  \item{\code{fit.rDoses.min}}{
#'    \link{numeric}: Lowest regenerative dose used for the fitting.}
#'  \item{\code{fit.rDoses.max}}{
#'    \link{numeric}: Highest regenerative dose used for the fitting.}
#' }
#' See also \link{calc_TL.MAAD.fit.Q} and \link{calc_TL.MAAD.fit.I}. \cr
#'
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
#' See also \link{plot_TL.MAAD}. \cr
#'
#' \bold{Rejection criteria} \cr
#' The rejection criteria are: \cr
#' \describe{
#'  \item{\code{testdose.error}}{
#'    \link{numeric}: Maximum error accepted on the testdose signal within the selected temperature interval.}
#'  \item{\code{paleodose.error}}{
#'    \link{numeric}: Maximum error accepted on the palaeodose signal within the selected temperature interval.}
#' }
#'
#' @return
#'  This function plots the results from the differents functions called using:
#'  \link{plot_extract.TL},
#'  \link{plot_remove.preheat}
#'  \link{plot_substract.background}
#'  \link{plot_align.peaks} and
#'  \link{plot_TL.MAAD}. \cr
#'
#'  This function saves a  file containing the luminescence data after the pretreatment in the specified folder.
#'
#'  Finally, it also provides an \link{list} containing: \cr
#'  \describe{
#'    \item{\code{De.GC}}{
#'      \link{data.frame}: Results obtained with the dose plateau approach and their uncertainties
#'      (\code{De}, \code{De.error}}
#'    \item{\code{De.DP}}{
#'      \link{data.frame}: Results obtained with the growth curve approach and their uncertainties
#'      (\code{De}, \code{De.error}}
#'  }
#'
#' @seealso
#'  \link{read_BIN2R},
#'  \link{Risoe.BINfileData2TLum.BIN.File},
#'  \link{mod_extract.TL},
#'  \link{mod_update.dType},
#'  \link{mod_remove.aliquot},
#'  \link{mod_remove.preheat},
#'  \link{mod_substract.background},
#'  \link{mod_align.peaks},
#'  \link{write_R2BIN},
#'  \link{TLum.BIN.File2TLum.Analysis} and
#'  \link{analyse_TL.MAAD}.
#'
#' @author David Strebler, University of Cologne (Germany), \cr David Strebler
#'
#' @export script_TL.MAAD

script_TL.MAAD <- function(

  file.name,

  eval.Tmin,

  eval.Tmax,

  relative.error= 0.05,

  remove.discs=NULL,

  file.parameters=list(file.extension =".binx",
                       folder.in = "./",
                       folder.out = "./"),

  aligning.parameters=list(peak.Tmin=NULL,
                           peak.Tmax=NULL,
                           no.testdose=FALSE),

  fitting.parameters=list(fit.method="LIN",
                          fit.weighted=FALSE,
                          fit.use.slope=FALSE,
                          fit.aDoses.min=0,
                          fit.aDoses.max=NA,
                          fit.rDoses.min=0,
                          fit.rDoses.max=NA),

  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA,
                           no.plot=FALSE),

  rejection.criteria = list(testdose.error = 10,
                            paleodose.error = 10)

){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if(missing(file.name)){
    stop("[script_TL.MAAD] Error: Input 'file.name' is missing.")

  }else if(!is.character(file.name)){
    stop("[script_TL.MAAD] Error: Input 'file.name' is not of type 'character'.")
  }

  if(missing(eval.Tmin)){
    stop("[script_TL.MAAD] Error: Input 'eval.Tmin' is missing.")
  }else if(!is.numeric(eval.Tmin)){
    stop("[script_TL.MAAD] Error: Input 'eval.Tmin' is not of type 'numeric'.")
  }

  if(missing(eval.Tmax)){
    stop("[script_TL.MAAD] Error: Input 'eval.Tmax' is missing.")
  }else if(!is.numeric(eval.Tmax)){
    stop("[script_TL.MAAD] Error: Input 'eval.Tmax' is not of type 'numeric'.")
  }

  if(!is.numeric(relative.error)){
    stop("[script_TL.MAAD] Error: Input 'relative.error' is not of type 'numeric'.")
  }

  if(!is.list(file.parameters)){
    stop("[script_TL.plateau] Error: Input 'file.parameters' is not of type 'list'.")
  }

  if(!is.list(aligning.parameters)){
    stop("[script_TL.MAAD] Error: Input 'aligning.parameters' is not of type 'list'.")
  }

  if(!is.list(fitting.parameters)){
    stop("[script_TL.MAAD] Error: Input 'fitting.parameters' is not of type 'list'.")
  }

  if(!is.list(plotting.parameters)){
    stop("[script_TL.MAAD] Error: Input 'plotting.parameters' is not of type 'list'.")
  }

  if(!is.list(rejection.criteria)){
    stop("[script_TL.MAAD] Error: Input 'rejection.criteria' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  # ------------------------------------------------------------------------------
  # Value check
  # ------------------------------------------------------------------------------

  # data recovery and pretreatment
  data <- script_TL.pretreatment(file.name = file.name,
                                 relative.error = relative.error,
                                 remove.discs = remove.discs,
                                 file.parameters = file.parameters,
                                 aligning.parameters = aligning.parameters,
                                 plotting.parameters = plotting.parameters)

  #Paleodose estimation

  # MAAD analysis
  temp.res <- analyse_TL.MAAD(object = data,
                              eval.Tmin=eval.Tmin,
                              eval.Tmax=eval.Tmax,
                              fitting.parameters=fitting.parameters,
                              plotting.parameters=plotting.parameters,
                              rejection.criteria=rejection.criteria
                              )

  temp.De.DP <- get_TLum.Results(object = temp.res, "DP")
  temp.De.GC <- get_TLum.Results(object = temp.res, "GC")

  De.values.DP <- data.frame(De=temp.De.DP$De,
                             De.error=temp.De.DP$De.error)
  De.values.GC <- data.frame(De=temp.De.GC$De,
                             De.error=temp.De.GC$De.error)

  result <- list(De.values.DP = De.values.DP,
                 De.values.GC = De.values.GC
  )

  return(result)
}
