#' The \code{LR.inter.cont} function performs a Wald test for an interaction term between a continuous variable and a SNP.
#' The function returns the p-value of the Wald test.
#' @export
#' @title Wald test for an interaction term
#' @name LR.inter.cont
#' @param x name or numeric vector corresponding to the SNP tested
#' @param formula an object of class "formula" : a symbolic description of the model to be fitted
#' without the interaction term
#' @param data a data frame containing the variables in the model
#' @param Z1 name of the variable which is tested in interaction with x (x:Z1)
#' @return pvalue of the Wald test for the interaction term
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}



LR.inter.cont <- function(x,formula,data,Z1){ 
  #   if(is.null(class.Z)) Z1 <- data[,Z] else Z1 <- factor(data[,Z])
  data <- data.frame(data,x=x,Z1=Z1)
  model1 <- glm(formula=update(formula,~.+x+x:Z1),data=data,family=binomial(link="logit"))
  pval <- summary(model1)$coef[dim(summary(model1)$coef)[1],4]
  return(pval) 
}








