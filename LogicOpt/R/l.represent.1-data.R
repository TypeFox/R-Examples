#' @title Logicopt truth table created from "represent.csv" dataset
#'
#' @description l.represent.1 is an logicopt compatible truth table generated 
#' from the QCA dataset "represent.csv" where output "WNP" is 1. 
#'
#' @docType data
#'
#' @usage data(l.represent.1)
#'
#' @format R data frame table 
#'
#' @keywords QCA datasets
#'
#' @source compass.org website and various QCA packages
#'
#' @examples
#' \dontrun{
#' # Read raw QCA dataset from csv file
#' inpath <- system.file("extdata/raw_qca/represent.csv", package="LogicOpt")
#' represent <- read.csv(inpath,row.names=1,na="")
#'
#' # Need to load a QCA package that contains truthTable function: (pick one)
#' # library(QCAGUI) 
#' # library(QCApro)
#'
#' # Create the QCA truth table 
#' q.represent.1 <- truthTable(represent, outcome = "WNP{1}")
#'
#' # Create the logicopt truth table
#' l.represent.1 <- QCAtt2LOtt(q.represent.1)
#' }
#'
#' # Load up truth table
#' data(l.represent.1)
#' 
#' # Optimize logicopt truth table  and print results
#' represent1 <- logicopt(l.represent.1,5,1,find_dc=TRUE,mode="multi-min")
#' print_multi_tt(represent1,eqn=TRUE,n_in=5,n_out=1,QCA=TRUE)
"l.represent.1"
