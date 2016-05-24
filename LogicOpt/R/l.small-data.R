#' @title Espresso truth table with 4 inputs and 3 outputs 
#'
#' @description Espresso compatible truth table generated 
#' from espresso format file small.esp.
#'
#' @docType data
#'
#' @usage data(l.small)
#'
#' @format R data frame table 
#'
#' @keywords Espresso truth-table 
#'
#' @examples
#' \dontrun{
#' # steps to recreate
#' inpath <- system.file("extdata/espresso/small.esp", package="LogicOpt")
#' l.small <- logicopt(esp_file=inpath,mode="echo")[1]
#' }
"l.small"
