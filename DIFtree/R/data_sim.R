#' @name data_sim
#' @title Simulated data set 
#' @description The data set is simulated from a Rasch model 
#' where some items exhibit differential item functioning. 
#' Existing differences in item difficulties are simulated by step-functions. 
#' The true, simulated DIF structure is described in Tutz and Berger (2015), Chapter 4.2. 
#' @docType data
#' @usage data(data_sim)
#' @format A data frame containing 500 observations on 5 variables:
#' 
#'  \bold{Y}   matrix with binary 0/1 response for 20 items
#'  
#'  \bold{x1}  binary covariate 1
#'
#'  \bold{x2}  metric covariate 1
#' 
#'  \bold{x3}  binary covariate 2
#' 
#'  \bold{x4}  metric covariate 2
#' 
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2015): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed Trees, Cornell University Library, arXiv:1511.07178
#' 
#' Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika, published online, DOI: 10.1007/s11336-015-9488-3 
#' 
#' @examples 
#' data(data_sim)
#'  
#' Y <- data_sim[,1]
#' X <- data_sim[,-1]
#' 
#' hist(rowSums(Y), breaks = 0:19 + 0.5)
#' summary(X)
#'   
#'  
NULL
