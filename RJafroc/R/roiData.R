#' @name roiData
#' @title An ROI dataset produced by a data simulator
#' @description This ROI dataset represents a simulation study in which two modalities are compared. 
#' There are 50 non-diseased and 40 diseased cases and each case is divided into 4 ROis. Five simulated 
#' readers rated each ROI in both modalities using a floating point scale, where larger ratings represented 
#' stronger confidence for presence of disease. The NL array contains the ratings of all non-diseased ROIs while
#' the LL array contains the ratings of all diseased ROIs. The *.xlsx file can be downloaded from
#' \url{http://www.devchakraborty.com/RoiData/roiData.xlsx}. The dataset file can then be read into a dataset 
#' object using \link{ReadDataFile} function.
#' @docType data
#' @usage roiData
#' @format
#' A dataset object consisting of a list containing 8 elements, see \link{RJafroc-package}.
#'
#' @source
#' The simulator used to generate the values is availabe from Dr. Chakraborty's website \url{http://www.devchakraborty.com/RoiData/RoiSimulator.zip}.
#' 
NULL 
