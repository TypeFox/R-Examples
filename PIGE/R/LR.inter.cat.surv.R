#' The \code{LR.inter.cat.surv} function performs a likelihood ratio test (LRT) for an interaction term between a categorical variable and a SNP (coded 0,1,2) in a cox model.
#' The function returns the p-value of the likelihood ratio test.
#' @export
#' @title Likelihhod ratio test for an interaction term
#' @name LR.inter.cat.surv
#' @param x name or numeric vector corresponding to the SNP tested.
#' @param formula an object of class "formula" : a symbolic description of the model to be fitted
#' without the interaction term.
#' @param data a data frame containing the variables in the model.
#' @param Z1 name of the variable which is tested in interaction with x (x:Z1).
#' @return p-value of the likelihood ratio test for the interaction term.
#' @author Benoit Liquet \email{benoit.liquet@@isped.u-bordeaux2.fr}\cr
#'  Therese Truong \email{therese.truong@@inserm.fr}


  LR.inter.cat.surv <- function(x,formula,data,Z1){ 
    if(is.character(x)) data <- data.frame(data,x=data[,x],Z1=Z1) else data <- data.frame(data,x=x,Z1=Z1)
    model1 <- coxph(formula=update(formula,~.+x+x:Z1),data=data)
    model2 <- coxph(formula=update(formula,~.+x),data=data)
    pval <- anova(model2,model1)[2,4]
    #   print(pval)
    return(pval)
    #return(model1)
  }