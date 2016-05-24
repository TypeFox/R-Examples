#' Plot empirical operating characteristic
#' 
#' Plot emprirical operating characteristics (operating points connected by straight lines) for specified trts and rdrs in the dataset.
#' 
#' @param dataset Dataset to be used for plotting.
#' @param trts List or vector: indices of modalities to be plotted. See "Details".
#' @param rdrs List or vector. indices of readers to be plotted. See "Details".
#' @param lgdPos The positioning of the legend: \code{"right"}(the default), \code{"left"}, \code{"top"} or \code{"bottom"}.
#' @param opChType Type of operating characteristic to be plotted. Available choices are \code{"ROC"}(the default), \code{"AFROC"} and \code{"FROC"}.
#' 
#' @details \strong{Note} that \code{trts} and \code{rdrs} are the vectors or list of \strong{indices} not \strong{IDs}. For example, if the ID of the 
#' first reader is "0". The corresponding value in \code{trts} should be \strong{1}  not \strong{0}. 
#' 
#' If both of \code{trts} and \code{rdrs} are vectors, all possible combinations will be plotted. 
#' 
#' If both of \code{trts} and \code{rdrs} are lists, they must have same length. Only the combination of modality and reader at same position will be 
#' plotted. If some elements of the lists are vectors, the averaged operating characteristic over them will be plotted. See "Examples".
#' 
#' @return A \pkg{ggplot2} object of the plotted operating characteristics and a data frame containing the points of the operating characteristics are returned. Following are the returned objects of 
#' \code{"ROC"} operating characteristics.
#' @return \item{ROCPlot}{\pkg{ggplot2} object: Use \code{print} function to display the saved object.}
#' @return \item{ROCPoints}{Data frame with four columns: abscissa, ordinate, class  (coding modality and reader) and type, 
#' which can be "individual" or "averaged".}
#' @examples
#' plotM <- c(1:2)
#' plotR <- c(1:3)
#' EmpiricalOpCharac(dataset = rocData, trts = plotM, rdrs = plotR,
#'                   lgdPos = "bottom", opChType = "ROC")
#' ## Above is the example of plotting individual ROC operating characteristics of modalities
#' ## 1 and 2 and readers 1 to 3. Six operating characteristics will be plotted, which are 
#' ## operating characteristics of reader 1 modality 1, reader 1 modality 2, reader 2 modality 
#' ## 1, reader 2 modality 2, reader 3 modality 1 and reader 3 modality 2.
#' 
#' plotM <- list(1, 2, c(1:2))
#' plotR <- list(2, c(2:3), c(1:3))
#' EmpiricalOpCharac(dataset = rocData, trts = plotM, rdrs = plotR,
#'                   lgdPos = "bottom", opChType = "ROC")
#' EmpiricalOpCharac(dataset = frocData, trts = plotM, rdrs = plotR,
#'                   lgdPos = "bottom", opChType = "AFROC")
#' EmpiricalOpCharac(dataset = frocData, trts = plotM, rdrs = plotR,
#'                   lgdPos = "bottom", opChType = "FROC")               
#' ## Above is the example of plotting three ROC, AFROC and FROC operating characteristics. 
#' ## They are the individual operating characteristic of modality 1 reader 2, the 
#' ## averaged operating characteristic of modality 2 and reader 2 and 3 and the averaged 
#' ## operating characteristic of modality 1 and 2 and reader 1 to 3.
#' 
#' @export
#' 
EmpiricalOpCharac <- function(dataset, trts, rdrs, lgdPos, opChType = "ROC") {
  if (dataset$dataType == "ROI")
    stop("The empirical operating characteristic of ROI dataset is unavailable.")
  if (opChType == "ROC") {
    return(PlotROC(dataset, trts, rdrs, lgdPos))
  } else if (opChType == "AFROC") {
    if (dataset$dataType == "ROC")
      stop("The AFROC operating characteristic of ROC dataset is unavailable.")
    return(PlotAFROC(dataset, trts, rdrs, lgdPos))
  } else if (opChType == "FROC") {
    if (dataset$dataType == "ROC")
      stop("The AFROC operating characteristic of ROC dataset is unavailable.")
    return(PlotFROC(dataset, trts, rdrs, lgdPos))
  } else {
    errMsg <- sprintf("%s is not an available operating characteristic type.", opChType)
    stop(errMsg)
  }
} 
