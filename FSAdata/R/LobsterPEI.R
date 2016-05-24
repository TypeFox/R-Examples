#' @title Catch and effort data for Prince Edward Island Lobster.
#' 
#' @description Catch and effort data for Lobster from 33 days in 1944 from the Tignish area of Prince Edward Island.
#' 
#' @details Catch (1000s of pounds) and effort (1000s of traps) of Lobster from 33 days in 1944 from the Tignish area of Prince Edward Island.  The data start on 2-May.  These data are from DeLury (1947) who used the data after 22-May (i.e., day 16) to illustrate his depletion method.  The data were also used in Example 7.1 of Seber (2002).  DeLury (1947) noted that the weight of Lobster did not change appreciably over time so that the poundage caught is a reasonable surrogate for numbers caught.
#' 
#' @name LobsterPEI
#' 
#' @docType data
#' 
#' @format A data frame with 34 observations on the following 3 variables.
#'  \describe{
#'    \item{day}{Day of the catch.  Day 1 is 2-May-1944.} 
#'    \item{catch}{Catch of Lobster in 1000s of pounds.} 
#'    \item{effort}{Total daily effort expended in 1000s of traps.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Population size 
#'    \item Abundance 
#'    \item Depletion methods
#'    \item Leslie method
#'    \item DeLury method 
#'    \item Catchability
#'  }
#'  
#' @concept Abundance 'Population Size' Leslie DeLury Depletion Catchability
#' 
#' @source From Table 1 of DeLury, D.B. 1947. On the estimation of biological populations. Biometrics 3:145-167.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(LobsterPEI)
#' str(LobsterPEI)
#' head(LobsterPEI)
#' 
NULL