#' @title Truth table from a Genetic Programming Use Case  
#'
#' @description This is an 8 input 3 output truth table from A.J. Keane 
#' that is used in a Genatic Programming use case for robot control.  The truth 
#' table has been modified slightly from the paper (referenced below) to specify 
#' Boolean outputs for the three possible output values of "zero", "one", 
#' and "minus" one.  This is slightly more informative than the example in the
#' paper.  
#'
#' @docType data
#'
#' @usage data(l.robot1)
#'
#' @format Espresso compatible truth table 
#'
#' @keywords Espresso truth-table Genetic Programming 
#'
#' @source
#' Keane, A.J. 2015. "Genetic Programming, Logic Design and Case-Based Reasoning for 
#' Obstacle Avoidance."  Learning and Intelligent Optimization: 9th International Conference, pages 104:118.
#'
#' @examples
#' \dontrun{
#' # steps to recreate l.robot1
#' inpath <- system.file("extdata/espresso/robot1_in.esp", package="LogicOpt")
#' l.robot1 <- logicopt(esp_file=inpath,mode="echo") 
#' }
#' 
#' # load l.robot1
#' data(l.robot1)
#'
#' # optimize l.robot1
#' robot1_opt <- logicopt(l.robot1,8,3)
#'
#' # optimized results have 13 rows that cover outputs zero, one, and minus
#' robot1_opt[2]
#'
#' # print optimized equations (where each output is 1)
#' print_multi_tt(robot1_opt,TRUE,8,3)

"l.robot1"
