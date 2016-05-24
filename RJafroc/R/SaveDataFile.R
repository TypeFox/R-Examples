#' Save ROC data file in a different format
#' 
#' Save ROC data file in a different format so it can be analyzed with alternate software.
#' 
#' @param dataset The dataset object to be saved in the specified format, see \link{RJafroc-package}.
#' @param fileName The file name of the output data file. The extension of the data file must match the corresponding format, see \link{RJafroc-package}.
#' @param format The format of the data file, which can be \code{"JAFROC"}, \code{"MRMC"} or \code{"iMRMC"}, see \link{RJafroc-package}.
#' @param dataDscrpt An optional string variable describing the data file, the default value is the variable name of \code{dataset}. 
#' The decription appears on the first line of *.lrc or *imrmc data file. This parameter is not used when saving dataset in other formats.
#' 
#' @examples
#' 
#' SaveDataFile(dataset = rocData, fileName = "rocData2.xlsx", format = "JAFROC")
#' SaveDataFile(dataset = rocData, fileName = "rocData2.csv", format = "MRMC")
#' SaveDataFile(dataset = rocData, fileName = "rocData2.lrc", format = "MRMC", 
#'              dataDscrpt = "ExampleROCdata")
#' SaveDataFile(dataset = rocData, fileName = "rocData2.txt", format = "MRMC", 
#'              dataDscrpt = "ExampleROCdata2")
#' SaveDataFile(dataset = rocData, fileName = "rocData.imrmc", format = "iMRMC", 
#'              dataDscrpt = "ExampleROCdata3") 
#' 
#' @importFrom tools file_ext
#' 
#' @export
#' 
SaveDataFile <- function(dataset, fileName, format = "JAFROC", dataDscrpt = paste0(deparse(substitute(dataset)), " Data File")) {
  if (format == "JAFROC") {
    return(SaveJAFROC(dataset, fileName))
  } else if (format == "iMRMC") {
    return(SaveImrmc(dataset, fileName, dataDscrpt))
  } else if (format == "MRMC") {
    if (file_ext(fileName) == "lrc") {
      return(SaveLrc(dataset, fileName, dataDscrpt))
    } else {
      return(SaveOrDbmMrmc(dataset, fileName))
    }
  } else {
    errMsg <- sprintf("%s is not an available file format.", format)
    stop(errMsg)
  }
} 
