#' Script for data pretreatment
#'
#' This script opens a .binx file and creates a \linkS4class{TLum.Analysis} object from it.
#' It just requires the name of the file with the TL curves and the relative error on the measurements.
#' It extracts the TL curves and updates the data types.
#'
#' @param file.name
#'  \link{character} (\bold{required}): Name of the file containing the luminescence data.
#' @param relative.error
#'  \link{numeric} (with default): Relative error of the TL signals.
#' @param protocol
#'  \link{character} (\bold{required}): Measurment protocol used.
#' @param file.parameters
#'  \link{list} (with default): list containing the file parameters. See details.
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
#' See also \link{plot_extract.TL}. \cr
#'
#' \bold{File parameters} \cr
#' The file parameters are:  \cr
#' \describe{
#'  \item{\code{file.extension}}{
#'    \link{character} (with default): extension of the file containing the luminescence data (.bin or .binx)}
#'  \item{\code{folder.in}}{
#'    \link{character} (with default): Folder containing the file with the luminescene data.}
#' }
#'
#' @return
#'  This function returns a \code{\linkS4class{TLum.Analysis}} object.
#'
#' @seealso
#'  \link{read_BIN2R},
#'  \link{Risoe.BINfileData2TLum.BIN.File},
#'  \link{TLum.BIN.File2TLum.Analysis},
#'  \link{mod_extract.TL},
#'  \link{mod_update.dType}.
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export script_TL.import

script_TL.import <- function(

  file.name,

  relative.error= 0.05,

  protocol = "Unknown",

  file.parameters=list(file.extension =".binx",
                       folder.in = "./"),

  plotting.parameters=list(plot.Tmin=0,
                           plot.Tmax=NA,
                           no.plot=FALSE)
){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if(missing(file.name)){
    stop("[script_TL.import] Error: Input 'file.name' is missing.")

  }else if(!is.character(file.name)){
    stop("[script_TL.import] Error: Input 'file.name' is not of type 'character'.")
  }

  if(!is.character(protocol)){
    stop("[script_TL.import] Error: Input 'protocol' is not of type 'character'.")
  }
  if(!is.numeric(relative.error)){
    stop("[script_TL.import] Error: Input 'relative.error' is not of type 'numeric'.")
  }

  if(!is.list(file.parameters)){
    stop("[script_TL.import] Error: Input 'file.parameters' is not of type 'list'.")
  }

  if(!is.list(plotting.parameters)){
    stop("[script_TL.import] Error: Input 'plotting.parameters' is not of type 'list'.")
  }
  # ------------------------------------------------------------------------------

  file.extension <- file.parameters$file.extension
  folder.in <- file.parameters$folder.in

  # ------------------------------------------------------------------------------
  # Value check
  if(relative.error > 1){
    warning("[script_TL.import] Warning: Input 'relative.error' > 1.")

  }else if(relative.error < -1){
    relative.error <- abs(relative.error)
    warning("[script_TL.import] Warning: Input 'relative.error' < -1.")

  }else if(relative.error < 0){
    relative.error <- abs(relative.error)
    warning("[script_TL.import] Warning: Input 'relative.error' < 0.")
  }

  if(!is.character(file.extension)){
    stop("[script_TL.import] Error: Input 'file.extension' is not of type 'character'.")
  }else if(file.extension !=  ".bin" &&  file.extension !=  ".binx"){
    stop("[script_TL.import] Error: Input 'file.extension' is not of '.bin' or '.binx'.")
    file.extension <- ".binx"
  }

  if(!is.character(folder.in)){
    warning("[script_TL.import] Error: Input 'folder.in' is not of type 'character'.")
    folder.in = "./"
  }
  # ------------------------------------------------------------------------------

  path.in <- paste(folder.in,file.name, file.extension,sep="")

  # Read file
  data.in <- read_BIN2R(path.in)

  data <- Risoe.BINfileData2TLum.BIN.File(object = data.in,
                                          relative.error = relative.error)

  data <- TLum.BIN.File2TLum.Analysis(object = data, protocol = protocol)

  #TL curve selection
  data <- mod_extract.TL(object = data, plotting.parameters = plotting.parameters)

  print("TL signals selected")

  #Identification of Preheat and Testdose
  data <- mod_update.dType(object = data)

  print("Preheat and tesdose identify")

  return(data)
}
