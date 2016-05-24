#' Script for the plateau test
#'
#' This script calls a series of data pretreatment functions before performing the plateau test.
#' It just requires the name of the file with the TL curves and the relative error on the measurements.
#'
#' @param file.name
#'  \link{character} (\bold{required}): Name of the file containing the luminescence data.
#' @param relative.error
#'  \link{numeric} (with default): Relative error of the TL signals.
#' @param remove.discs
#'  \link{numeric}  (with default): list containing the position of the aliquots that shall be removed.
#' @param file.parameters
#'  \link{list} (with default): list containing the file parameters. See details.
#' @param aligning.parameters
#'  \link{list} (with default): list containing the aligning parameters. See details.
#' @param plotting.parameters
#'  \link{list} (with default): list containing the plotting parameters. See details.
#'
#' @details
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
#'
#' @return
#'  This function plots the results from the differents functions called using:
#'  \link{plot_extract.TL},
#'  \link{plot_remove.preheat}
#'  \link{plot_substract.background}
#'  \link{plot_align.peaks} and
#'  \link{plot_TL.plateau}. \cr
#'
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
#'  \link{analyse_TL.plateau}.
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export script_TL.plateau

script_TL.plateau <- function(

  file.name,


  relative.error= 0.05,

  remove.discs=NULL,

  file.parameters=list(file.extension =".binx",
                       folder.in = "./",
                       folder.out = "./"),

  aligning.parameters=list(peak.Tmin=NULL,
                           peak.Tmax=NULL,
                           no.testdose=FALSE),


  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA,
                           plateau.Tmin=0,
                           plateau.Tmax=0,
                           no.plot=FALSE)


){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if(missing(file.name)){
    stop("[script_TL.plateau] Error: Input 'file.name' is missing.")
  }else if(!is.character(file.name)){
    stop("[script_TL.plateau] Error: Input 'file.name' is not of type 'character'.")
  }

  if(!is.numeric(relative.error)){
    stop("[script_TL.plateau] Error: Input 'relative.error' is not of type 'numeric'.")
  }

  if(!is.list(file.parameters)){
    stop("[script_TL.plateau] Error: Input 'file.parameters' is not of type 'list'.")
  }

  if(!is.list(aligning.parameters)){
    stop("[script_TL.plateau] Error: Input 'aligning.parameters' is not of type 'list'.")
  }

  if(!is.list(plotting.parameters)){
    stop("[script_TL.plateau] Error: Input 'plotting.parameters' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  data <- script_TL.pretreatment(file.name = file.name,
                                 relative.error = relative.error,
                                 remove.discs = remove.discs,
                                 file.parameters = file.parameters,
                                 aligning.parameters = aligning.parameters,
                                 plotting.parameters = plotting.parameters)

  # MAAD analysis
  temp.res <- analyse_TL.plateau(object = data,
                                 plotting.parameters=plotting.parameters)

  return(temp.res)
}
