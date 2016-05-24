#' OR and their confidence intervals for logistic regressions
#'
#' Computes odd ratios and their confidence intervals for logistic regressions
#'
#' @param model a \code{glm} object
#' @param alpha type I error, 0.05 by default
#' @return A matrix with the estimaed coefficients of the logistic model, their s.e., z-values, p-values, OR and CI of the OR
#' @author Hugo Varet
#' @examples
#' IC_OR_glm(glm(inherit~sex+age,data=cgd,family="binomial"))

IC_OR_glm=function(model,alpha=0.05){
  # model doit etre un objet résultant de glm
  tab=matrix(nrow=nrow(summary(model)$coefficients),ncol=ncol(summary(model)$coefficients)+3,
           dimnames = list(c(rownames(summary(model)$coefficients)),c(colnames(summary(model)$coefficients),"OR","IC.inf","IC.sup")))
  tab[,1:ncol(summary(model)$coefficients)]=round(summary(model)$coefficients,digits=3)
  for (i in 1:(nrow(summary(model)$coefficients))){
    # i = indice de la covariable
    OR=round(exp(summary(model)$coefficients[i,1]),digits=3)
    ICinf=round(exp(summary(model)$coefficients[i,1]-qnorm(1-alpha/2)*summary(model)$coefficients[i,2]),digits=3)
    ICsup=round(exp(summary(model)$coefficients[i,1]+qnorm(1-alpha/2)*summary(model)$coefficients[i,2]),digits=3)
    tab[i,5]=as.numeric(OR)
    tab[i,6]=ICinf
    tab[i,7]=ICsup
  }
  p=tab[,4]
  signif=ifelse(p>=0.1,"",ifelse(p>=0.05,".",ifelse(p>=0.01,"*",ifelse(p>=0.001,"**","***"))))
  res=cbind(tab,signif)
  colnames(res)=c(colnames(tab),"")
  return(noquote(res))
}

#y=rbinom(1000,1,0.4)
#x1=runif(1000)
#x2=runif(1000)
#x3=runif(1000)
#x4=runif(1000)
#
#model=glm(y~x1,family="binomial")
#IC_OR_glm(model)
