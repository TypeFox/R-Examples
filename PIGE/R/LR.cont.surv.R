#' The \code{LR.cont.surv} function performs a Wald test for a SNP effect based on a cox model.
#' The function returns the p-value of the Wald test. It is an Internal function used by the \link{permutation.snp} function.
#' @export
#' @title Wald test for an adjusted model.
#' @name LR.cont.surv
#' @param x name or numeric vector corresponding to the SNP tested.
#' @param formula an object of class "formula" : a symbolic description of the model to be fitted
#' without the interaction term.
#' @param data a data frame containing the variables in the model.
#' @return p-value of the Wald test for a SNP effect.
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}

LR.cont.surv <- function(x,formula,data){ 
  data <- data.frame(data,x=x)
  model1 <- coxph(formula=update(formula,~.+x),data=data)
  pval <- summary(model1)$coef[dim(summary(model1)$coef)[1],5]
  return(pval) 
}



