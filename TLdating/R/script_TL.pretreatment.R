#' Script for data pretreatment
#'
#' This script call a series of data pretreatment functions for TL dating.
#' It only requires the name of the files with the TL curves and the relative error on the measurement.
#'
#' @param file.name
#'  \link{character} (\bold{required}): Name of the file containing the luminescence data.
#' @param relative.error
#'  \link{numeric} (with default): Relative error of the TL signals.
#' @param remove.discs
#'  \link{numeric}  (with default): list containing the position of the aliquots to remove.
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
#'    \link{numeric}: Lower boundary for looking at the peak maximum position.}
#'  \item{\code{peak.Tmax}}{
#'    \link{numeric}: Upper boundary for looking at the peak maximum position.}
#'  \item{\code{no.testdose}}{
#'    \link{logical}: If \code{TRUE}, the function will use the Lx curves rather the Tx curves as reference for the peak maximum position.}
#'  }
#'
#' \bold{Plotting parameters} \cr
#' The plotting parameters are:  \cr
#' \describe{
#'  \item{\code{plot.Tmin}}{
#'    \link{numeric}: Lower temperature plotted.}
#'  \item{\code{plot.Tmax}}{
#'    \link{numeric}: Higher temperature plotted.}
#'  \item{\code{no.plot}}{
#'    \link{logical}: If \code{TRUE}, the results will not be plotted.}
#' }
#' See also \link{plot_TL.MAAD}. \cr
#'
#' \bold{File parameters} \cr
#' The plotting parameters are:  \cr
#' \describe{
#'  \item{\code{file.extension}}{
#'    \link{character} (with default): extension of the file containing the luminescence data (.bin or .binx)}
#'  \item{\code{folder.in}}{
#'    \link{character} (with default): Folder containing the file with the luminescene data.}
#'  \item{\code{folder.out}}{
#'    \link{character} (with default): Folder containing the file with the new luminescene data.}
#' }
#' see also \link{mod_update.dType}.
#'
#' @return
#'  This function return a \code{\linkS4class{TLum.Analysis}} where the preheat were removed, the background substract and the peaks aligned.
#'  Its  save the result as a .binx file il the specified folder.
#'  And, its plots the results from the differents functions called using:
#'  \link{plot_extract.TL},
#'  \link{plot_remove.preheat},
#'  \link{plot_substract.background} and
#'  \link{plot_align.peaks}. \cr
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
#'  \link{write_R2BIN}.
#'
#' @author David Strebler, University of Cologne (Germany), \cr David Strebler
#'
#' @export script_TL.pretreatment

script_TL.pretreatment <- function(

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
                           no.plot=FALSE)

){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if(missing(file.name)){
    stop("[script_TL.pretreatment] Error: Input 'file.name' is missing.")

  }else if(!is.character(file.name)){
    stop("[script_TL.pretreatment] Error: Input 'file.name' is not of type 'character'.")
  }

  if(!is.numeric(relative.error)){
    stop("[script_TL.pretreatment] Error: Input 'relative.error' is not of type 'numeric'.")
  }

  if(!is.list(file.parameters)){
    stop("[script_TL.pretreatment] Error: Input 'plotting.parameters' is not of type 'list'.")
  }

  if(!is.list(aligning.parameters)){
    stop("[script_TL.pretreatment] Error: Input 'aligning.parameters' is not of type 'list'.")
  }

  if(!is.list(plotting.parameters)){
    stop("[script_TL.pretreatment] Error: Input 'plotting.parameters' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  folder.out <- file.parameters$folder.out
  file.extension <- file.parameters$file.extension

  # ------------------------------------------------------------------------------
  # Check Value
  if(!is.character(folder.out)){
    warning("[script_TL.pretreatment] Error: Input 'folder.out' is not of type 'character'.")
    folder.out = "./"
  }

  if(!is.character(file.extension)){
    stop("[script_TL.pretreatment] Error: Input 'file.extension' is not of type 'character'.")
  }else if(file.extension !=  ".bin" &&  file.extension !=  ".binx"){
    stop("[script_TL.pretreatment] Error: Input 'file.extension' is not of '.bin' or '.binx'.")
    file.extension <- ".binx"
  }

  # ------------------------------------------------------------------------------

  # TL curve recovery
  data <- script_TL.import(file.name = file.name,
                          relative.error = relative.error,
                          file.parameters = file.parameters,
                          plotting.parameters = plotting.parameters)

  #Problematic aliquots removal
  if(!is.null(remove.discs)){
    data <- mod_remove.aliquot(object = data,
                               list = remove.discs)
    print(paste("Aliquot", remove.discs, "removed"))
  }

  # Preheat removal
  data <- mod_remove.preheat(object = data,
                             plotting.parameters = plotting.parameters)
  print("Preheat removed")

  # Background substraction
  data <- mod_substract.background(object = data)
  print("Background substracted")

  # Peaks alignement
  data <- mod_align.peaks(object=data,
                          aligning.parameters=aligning.parameters,
                          plotting.parameters=plotting.parameters)
  print("Peaks Shifted")

  #Saving of preliminary results

  file.out <- paste("new_",file.name,sep="")
  script_TL.export(object = data,
                   file.name = file.out,
                   file.parameters = file.parameters)

  # path.out <-  paste(folder.out,"new_",file.name,file.extension,sep="")
  #
  # data.out <- TLum.Analysis2TLum.BIN.File(data)
  # data.out <- TLum.BIN.File2Risoe.BINfileData(data.out)
  #
  # write_R2BIN(object = data.out,
  #             file =  path.out)

  print("File saved")

  return(data)
}
