#' Reads the data file and creates a dataset object
#' 
#' Reads the dataset file to be analyzed and creates a dataset object for subsequent analysis.
#' 
#' @param fileName A string specifying the name of the file that contains the dataset. 
#' The extension of the file must match the corresponding format specified below.
#' @param format A string specifying the format of the data file. It can be \code{"JAFROC"} (the default), \code{"MRMC"} or \code{"iMRMC"}. 
#' For \code{"MRMC"} the format is determined by the extension of the data file as specified in \url{http://perception.radiology.uiowa.edu/}: 
#' the file extension can be .csv or .txt or .lrc. For file extension .imrmc the format is described in \url{https://code.google.com/p/imrmc/}.
#' @param delimiter The string delimiter to be used for the \code{"MRMC"} format ("," is the default), see \url{http://perception.radiology.uiowa.edu/}.
#' This parameter is not used when reading \code{"JAFROC"} or \code{"iMRMC"} data files.
#' 
#' @return A dataset with the specified structure, see \link{RJafroc-package}.
#' 
#' @examples
#' 
#' fileName <- system.file("tests", "rocData.xlsx", package = "RJafroc")
#' RocDataXlsx<- ReadDataFile(fileName)
#' 
#' \dontrun{
#' fileName <- system.file("tests", "rocData.csv", package = "RJafroc")
#' RocDataCsv<- ReadDataFile(fileName, format = "MRMC")
#' 
#' fileName <- system.file("tests", "rocData.imrmc", package = "RJafroc")
#' RocDataImrmc<- ReadDataFile(fileName, format = "iMRMC")
#' 
#' fileName <- system.file("tests", "frocData.xlsx", package = "RJafroc")
#' FrocDataXlsx <- ReadDataFile(fileName)
#' 
#' fileName <- system.file("tests", "roiData.xlsx", package = "RJafroc")
#' RoiDataXlsx <- ReadDataFile(fileName)
#' }
#' 
#' @importFrom tools file_ext
#' 
#' @export
#' 
ReadDataFile <- function(fileName, format = "JAFROC", delimiter = ",") {
  if (format == "JAFROC") {
    if (!(file_ext(fileName) %in% c("xls", "xlsx"))) 
      stop("The extension of JAFROC data file must be \"*.xls\" or \"*.xlsx\" ")
    return(ReadJAFROC(fileName))
  } else if (format == "iMRMC") {
    return(ReadImrmc(fileName))
  } else if (format == "MRMC") {
    if (file_ext(fileName) == "lrc") {
      return(ReadLrc(fileName))
    } else {
      return(ReadOrDbmMrmc(fileName, delimiter))
    }
  } else {
    errMsg <- sprintf("%s is not an available file format.", format)
    stop(errMsg)
  }
} 
