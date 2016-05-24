#' Problem: Comparison between several Machine Learning algorithms from the Weka library
#'
#' A problem containing experimental data obtaining by comparing several instances
#' of Machine Algorithms from the Weka library. The variables are as follows:
#'
#' \itemize{
#'   \item method. Classification algorithms used in the experimen (NaiveBayes, J48, IBk)
#'   \item problem. Problems used as benchmark in the comparison, up to 12.
#'   \item featureSelection. Boolean parameter indicating if the data was preprocessed
#'   \item fold. For each configuration a 10-fold cross validation was performed. This variable is a numeric value ranging from 1 to 10.
#'   \item accuracy. This is a measure of the performance of each algorithm. Representing the percentage of correctly classified instances.
#'   \item trainingTime. A second measure of performance. This one indicates the time in seconds that took the algorithm to build the model.
#' }
#' @docType data
#' @keywords problems
#' @name wekaExperiment
#' @usage data(wekaExperiment)
#' @format A data frame with 
#'
"wekaExperiment"
