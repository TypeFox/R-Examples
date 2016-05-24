#' @name frocData
#' @title An example FROC dataset provided by Dr. Federica Zanca.
#' @description In this example dataset there are 100 non-diseased and 100 lesion containing mammography images. The maximum number
#' of lesions per diseased case is 3. Four radiologists intepreted all images using five image processing algorithms (modalities). 
#' They marked and rated suspicious regions on an integer 1 - 5 scale, where larger ratings represented greater 
#' confidence in presence of disease. Only modalities 4 and 5 have been retained in this dataset. The *.xlsx file can be 
#' downloaded from \url{http://www.devchakraborty.com/FrocData/frocData.xlsx}. The dataset file can then be read into a dataset 
#' object using \link{ReadDataFile} function.
#' @docType data
#' @usage frocData
#' @format
#' A dataset object consisting of a list containing 8 elements, see \link{RJafroc-package}.
#' 
#' @source
#' Zanca, F., Chakraborty, D. P., Van Ongeval, C., Jacobs, J., Claus, F., Marchal, G., & Bosmans, H. (2008). 
#' An improved method for simulating microcalcifications in digital mammograms. Medical Physics, 35(9), 4012-8.
#' 
NULL 
