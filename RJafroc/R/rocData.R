#' @name rocData
#' @title An ROC dataset originally provided by Dr. Kevin Berbaum, U of Iowa, ca. 2002.
#' @description This ROC dataset was collected by Carolyn Van Dyke, MD, in an ROC study that compared the performance of two modalities 
#' (Spin Echo MRI and cine MRI). There are 69 non-diseased and 45 diseased cases. Five radiologists interpreted all images in both modalities
#' using an integer 1 - 5 ratings scale, where larger ratings represented stronger confidence for presence of disease. 
#' The *.xlsx file can be downloaded from \url{http://www.devchakraborty.com/RocData/rocData.xlsx}. The dataset file can then 
#' be read into a dataset object using \link{ReadDataFile} function. The *.csv and *.lrc files can be downloaded from OR-DBM MRMC 
#' website, currently \url{http://perception.radiology.uiowa.edu/}. Alternatively they can be downloaded from 
#' \url{http://www.devchakraborty.com/RocData/rocData.csv} or \url{http://www.devchakraborty.com/RocData/rocData.lrc}. The iMRMC file 
#' can be downloaded from \url{http://www.devchakraborty.com/RocData/rocData.imrmc}.  
#' @docType data
#' @usage rocData
#' @format
#' A dataset object consisting of a list containing 8 elements, see \link{RJafroc-package}.
#'
#' @source
#' Van Dyke, C. W., et al. "Cine MRI in the diagnosis of thoracic aortic dissection." 79th Annual Meeting of the Radiological 
#' Society of North America, Radiological Society of North America, Chicago, Illinois. 1993.
#' 
NULL 
