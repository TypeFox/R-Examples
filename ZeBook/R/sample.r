################################################################################
# "Working with dynamic models for agriculture"
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title Generate a random plan as a data frame. Columns are parameters. Values have uniform distribution
#' @description according to minimal and maximal values defined in a model.factors matrix 
#' @param model.factors : matrix defining minimal (binf) and maximal values (bsup) for a set of p parameters
#' @param N : size of sample
#' @return parameter matrix of dim = (N, p)
#' @export
param.runif=function(model.factors,N){ 
  X=data.frame("id"=1:N)
  for (p in colnames(model.factors)) {
      X=cbind(X,runif(N,min = model.factors["binf",p], max = model.factors["bsup",p]))
  } 
  X=X[,-1]
  names(X)=colnames(model.factors)
  return(X)
}
# end of function
################################################################################
#' @title Generate a random plan as a data frame. Columns are parameters. Values have triangle distribution
#' @description according to nominal, minimal and maximal values defined in a model.factors matrix
#' @param model.factors : matrix defining nominal, minimal (binf), maximal values (bsup) for a set of p parameters
#' @param N : size of sample
#' @return parameter matrix of dim = (N, p)
#' @export
param.rtriangle = function(model.factors,N)
{
  X=data.frame("id"=1:N)
  for (p in colnames(model.factors)) {
      # rtriangle : a : lower limit, b : upper limit and c : mode of the distribution
      X=cbind(X,rtriangle(N,a = model.factors["binf",p], b = model.factors["bsup",p], c = model.factors["nominal",p]))
  }
  X=X[,-1]
  names(X)=colnames(model.factors)
  return(X)
}
# end of function
################################################################################
#' @title Build the q.arg argument for the  FAST function (sensitivity analysis)
#' @description according to minimal and maximal values defined in a model.factors matrix 
#' @param model.factors : matrix defining minimal (binf) and maximal values (bsup) for a set of p parameters
#' @return a list of list
#' @export 
q.arg.fast.runif=function(model.factors){
  q.arg.fast = list()
  for (p in colnames(model.factors)) {
    q.arg.fast=c(q.arg.fast,list(list(min = model.factors["binf",p], max = model.factors["bsup",p])))
  }
  return(q.arg.fast) 
}
# end of function
