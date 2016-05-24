################################################################
### utility functions for the sga project
#' @name utility
#' @rdname utility
#' @aliases utility
#' @title Utility measures 
#' @description Various utility measues that basically compares two data sets
#' @param x a data.frame, typically the original data set. For \code{utilityIndicator} this should be a vector of length 1.
#' @param y a data.frame, typically the corresponding synthetic data set. For \code{utilityIndicator} this should be a vector of length 1.
#' @param type which measure
#' \itemize{
#'  \item{compareColumns}{compares the intersection of variables}
#'  \item{compareRows}{compares the number of rows}
#'  \item{compareRowsHH}{compares the number of housholds}
#'  \item{compareNA}{compares the number of missings}
#' }
#' @author Matthias Templ, Maxime Bergeaut
NULL

#' @describeIn utility comparisons of two data sets
#' @param hhid index or name of variable containing the houshold ID
#' @return the measure(s) of interest
#' @examples 
#' data(eusilcS)
#' data(eusilcP)
#' ## for fast caluclations, took a subsample
#' eusilcP <- eusilcP[1:15000, ]
#' utility(eusilcS, eusilcP)
#' @export
# Basic similarity measures about structure of PUFs
utility <- function(x, y, type=c("all", "compareColumns", "compareRows", "compareRowsHH", "compareNA"), hhid = NULL){
  type <- match.arg(type)
  if(type=="all" | type=="compareColumns"){
    # Suppression of variables with missing values for every record
    indNAx <- apply(x, 2, function(x) all(is.na(x)))
    indNAy <- apply(y, 2, function(x) all(is.na(x))) 
    colnames(x)[!indNAx] %in% colnames(x)[!indNAy]
    common <- intersect(colnames(x)[!indNAx], colnames(y)[!indNAy])
    compareColumns <- length(common) / length(indNAx)
  } else {
    compareColumns <- NULL
  }
  
  if(type=="all" | type=="compareRows"){
    compareRows <- nrow(x) / nrow(y)
  } else {
    compareRows <- NULL
  }
  
  if(type=="all" | type=="compareRowsHH"){
    DT <- data.table(x)
    setkeyv(DT, hhid)
    DTy <- data.table(y)
    setkeyv(DTy, hhid)    
    compareRowsHH <- nrow(unique(DT)) / nrow(unique(DTy))
  } else {
    compareRowsHH <- NULL
  }
  
  if(type=="all" | type=="compareNA"){
    # Comparison for the NA structure for variables present in both PUF and SUF.
    listCommonVariables <- intersect(colnames(x),colnames(y))
    x_SimilarVariables <- subset(x, select=listCommonVariables)
    y_SimilarVariables <- subset(y, select=listCommonVariables)
    
    puf <- sum(is.na(x_SimilarVariables))
    suf <- sum(is.na(y_SimilarVariables))
    
    if(suf > 0 & puf >= suf ){ 
      compareNA <- suf / puf
    } 
    else if(suf > 0 & puf < suf) {
      compareNA <- puf / suf  
    }
    else if(suf == 0 & puf == 0) {
      compareNA <- 1  
    } 
    else if(suf > 0 & puf == 0){
      compareNA <- 0
    } 
  } else {
    compareNA <- NULL
  }
  comparisons <- list("compareColumns"=compareColumns,
                   "compareRows"=compareRows,
                   "compareRowsHH"=compareRowsHH,
                   "compareNA"=compareNA)
  return(comparisons)
}

#' @describeIn utility comparison of number of categories
#' @param varx name or index of a variable in data.frame x
#' @param vary NULL or name or index of a variable in data.frame y corresponding to variable varx in data.frame x. 
#' If NULL, the names of the selected variable should be the same in both x and y.
#' @examples 
#' data(eusilcS)
#' data(eusilcP)
#' utilityModal(eusilcS, eusilcP, "age")
#' utilityModal(eusilcS, eusilcP, "pl030", "ecoStat")
#' @export
utilityModal <- function(x, y, varx, vary=NULL){
  if(is.null(vary)){ 
    vary <- varx
  }
  measure5 <- length(table(x[, varx])) / length(table(y[, vary]))
  return("utiltiyModal" = measure5)
}

#' @describeIn utility difference between two values
#' @examples 
#' data(eusilcS)
#' data(eusilcP)
#' m1 <- meanWt(eusilcS$age, eusilcS$rb050) 
#' m2 <- mean(eusilcP$age)
#' utilityIndicator(m1, m2)
#' @export
utilityIndicator <- function(x, y){
  measure6 <- abs(x - y) / y 
  return("utilityIndicator"=measure6)
}

