################################################################################
# "Working with dynamic models for agriculture"
# R script for practical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013/06/14
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################################################################
#' @title The SeedWeight model
#' @description The SeedWeight model is a logistic model of grain weight over time in wheat. The model was proposed by Darroch & Baker (1990) in a study of grain filling in three spring wheat genotypes.
#'  This model has a single input variable, degree days after anthesis noted DD, and three parameters, noted W, B and C.
#' Parameters are estimated from observations.
#' @param DD : degree days after anthesis
#' @param W : parameter of the model
#' @param B : parameter of the model
#' @param C : parameter of the model
#' @return Seed Weight for each TT
#' @examples plot(1:500,seedweight.model(1:500, W=30,B=4,C=0.020),type="l",
#'    xlab="degree days after anthesis", ylab="grain weight")
#' @export
seedweight.model<-function(DD, W, B, C){
return(W/(1 + exp(B - C * DD)))}

# End of file
