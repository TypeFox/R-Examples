#' Imputes missing standard deviations in a dataset.
#'
#' Imputes (fills gaps) of missing standard deviations (SD) using simple imputation
#' methods following Bracken (1992) and Rubin and Schenker's (1991) "hot deck" 
#' approach.
#'
#' @param aDataFrame A data frame containing columns with missing SD's (coded as 
#' \code{NA}) and their complete means (used only for nearest-neighbor method).
#' @param columnSDnames Label of the column(s) with missing SD.  Can be a string 
#'    or list of strings.  
#' @param columnXnames Label of the column(s) with means (X) for each SD.  Can be 
#'    a string or list of strings.  Must be complete with no missing data.
#' @param method The method used to impute the missing SD's.  The default is 
#'    \code{"Bracken1992"} which applies Bracken's (1992) approach to impute SD using
#'    the coefficient of variation from all complete cases.  Other options include:
#'    \code{"HotDeck"} which applies Rubin and Schenker's (1991) resampling approach to
#'    fill gaps of missing SD from the SD's with complete information, and 
#'    \code{"HotDeck_NN"} which resamples from complete cases with means that are similar
#'    to missing SD's.
#'  @param range A positive number on the range of neighbours to sample from for
#'    imputing SD's.  Used in combination with \code{"HotDeck_NN"}. The default 
#'    is 3; which indicates that the 3 means that are most similar in rank order
#'    to the mean with the missing SD will be resampled. 
#'  @param M The number of imputed datasets to return.  Currently only works
#'    for \code{"HotDeck"} method.
#'
#' @return An imputed (complete) dataset. 
#'
#' @references Bracken, M.B. 1992. Statistical methods for analysis of effects
#'    of treatment in overviews of randomized trials. Effective care of the
#'    newborn infant (eds J.C. Sinclair and M.B. Bracken), pp. 
#'    13-20. Oxford University Press, Oxford.
#' @references Rubin, D.B. and Schenker, N. 1991. Multiple imputation in 
#'    health-care databases: an overview and some applications.  Statistics 
#'    in Medicine 10: 585-598.
#'
#' @export impute_SD

impute_SD <- function(aDataFrame, 
                      columnSDnames, 
                      columnXnames, 
                      method = "Bracken1992",
                      range = 3, 
                      M = 1) {

  if(method == "Bracken1992") {
    imputedData <- impute_SD_Bracken1992(aDataFrame,
                                         columnSDnames,
                                         columnXnames)
  } else if (method == "HotDeck") {
    imputedData <- impute_SD_HotDeck_fullRandom(aDataFrame,
                                                columnSDnames,
                                                M)
  } else { 
    imputedData <- impute_SD_HotDeck_nearestNeighbour(aDataFrame,
                                                      columnSDnames,
                                                      columnXnames,
                                                      range)
  }
                 
  return(imputedData)
}

