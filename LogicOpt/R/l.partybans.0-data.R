#' @title Logicopt truth table created from "partybans.csv" dataset
#'
#' @description l.partybans.0 is an logicopt compatible truth table generated 
#' from the QCA dataset "partybans.csv" where output "PB" is 0. 
#'
#' @docType data
#'
#' @usage data(l.partybans.0)
#'
#' @format R data frame table 
#'
#' @keywords QCA datasets
#'
#' @source compass.org website
#'
#' @examples
#' \dontrun{
#' # Read raw QCA dataset from csv file
#' inpath <- system.file("extdata/raw_qca/partybans.csv", package="LogicOpt")
#' partybans <- read.csv(inpath,row.names=1,na="")
#'
#' # Load QCA package 
#' library(QCAGUI) 
#'
#' # Create the QCA truth table 
#' q.partybans.0 <- truthTable(partybans, conditions = c("C","F","T","R","V"), outcome = "PB{0}")
#'
#' # Create the logicopt  truth table
#' l.partybans.0 <- QCAtt2LOtt(q.partybans.0)
#' }
#'
#' # Load up logicopt truth table 
#' data(l.partybans.0)
#' 
#' # Optimize and print logicopt truth table
#' partybans0 <- logicopt(l.partybans.0,5,1,find_dc=TRUE,mode="multi-min")
#' print_multi_tt(partybans0,eqn=TRUE,n_in=5,n_out=1,QCA=TRUE)
"l.partybans.0"
