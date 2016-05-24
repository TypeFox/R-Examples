#' The \code{LR.cont} function performs a Wald test for a SNP effect.
#' The function returns the p-value of the Wald test. It is an Internal function used by the \link{permutation.snp} function.
#' @export
#' @title Wald test for an adjusted model.
#' @name LR.cont
#' @param x name or numeric vector corresponding to the SNP tested.
#' @param formula an object of class "formula" : a symbolic description of the model to be fitted
#' without the interaction term.
#' @param data a data frame containing the variables in the model.
#' @return p-value of the Wald test for a SNP effect.
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}

LR.cont <- function(x,formula,data){ 
  data <- data.frame(data,x=x)
  model1 <- glm(formula=update(formula,~.+x),data=data,family=binomial(link="logit"))
  pval <- summary(model1)$coef[dim(summary(model1)$coef)[1],4]
  return(pval) 
}



