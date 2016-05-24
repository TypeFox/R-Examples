#' Script for data export
#'
#' This script creates a .binx file from a \linkS4class{TLum.Analysis} object.
#' It just requires the name of the file and the \linkS4class{TLum.Analysis} object.
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the TL curves to export.
#' @param file.name
#'  \link{character} (\bold{required}): Name of the file containing the luminescence data.
#' @param file.parameters
#'  \link{list} (with default): list containing the file parameters. See details.
#'
#'
#' \bold{File parameters} \cr
#' The file parameters are:  \cr
#' \describe{
#'  \item{\code{file.extension}}{
#'    \link{character} (with default): extension of the file containing the luminescence data (.bin or .binx)}
#'  \item{\code{folder.out}}{
#'    \link{character} (with default): Folder containing the file with the luminescene data.}
#' }
#'
#' @return
#'  This function returns a \code{\linkS4class{TLum.Analysis}} object.
#'
#' @seealso
#'  \link{write_R2BIN},
#'  \link{TLum.BIN.File2Risoe.BINfileData},
#'  \link{TLum.Analysis2TLum.BIN.File}.
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export script_TL.export

script_TL.export <- function(

  object,

  file.name,

  file.parameters=list(file.extension =".binx",
                       folder.out = "./")
){
  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if(missing(file.name)){
    stop("[script_TL.import] Error: Input 'file.name' is missing.")

  }else if(!is.character(file.name)){
    stop("[script_TL.import] Error: Input 'file.name' is not of type 'character'.")
  }

  if(!is.list(file.parameters)){
    stop("[script_TL.import] Error: Input 'file.parameters' is not of type 'list'.")
  }

  # ------------------------------------------------------------------------------

  file.extension <- file.parameters$file.extension
  folder.out <- file.parameters$folder.out

  # ------------------------------------------------------------------------------
  # Value check

  if(!is.character(file.extension)){
    stop("[script_TL.import] Error: Input 'file.extension' is not of type 'character'.")
  }else if(file.extension !=  ".bin" &&  file.extension !=  ".binx"){
    stop("[script_TL.import] Error: Input 'file.extension' is not of '.bin' or '.binx'.")
    file.extension <- ".binx"
  }

  if(!is.character(folder.out)){
    warning("[script_TL.import] Error: Input 'folder.out' is not of type 'character'.")
    folder.out = "./"
  }
  # ------------------------------------------------------------------------------

  path.out <- paste(folder.out,file.name, file.extension,sep="")


  # write file

  data <- TLum.Analysis2TLum.BIN.File(object)

  data.out <- TLum.BIN.File2Risoe.BINfileData(data)

  write_R2BIN(object = data.out, file =  path.out)
}
