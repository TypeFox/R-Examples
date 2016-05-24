#' RR and their confidence intervals for Cox models
#'
#' Computess risk ratios and their confidence intervals for Cox models
#'
#' @param model a \code{coxph} object
#' @param alpha type I error, 0.05 by default
#' @param sided 1 or 2 for one or two-sided
#' @return A matrix with the estimaed coefficients of the Cox model, their s.e., z-values, p-values, RR and CI of the RR
#' @author Hugo Varet
#' @examples
#' cgd$time=cgd$tstop-cgd$tstart
#' IC_RR_coxph(coxph(Surv(time,status)~sex+age,data=cgd),alpha=0.05,sided=1)

IC_RR_coxph=function(model,alpha=0.05,sided=2){
  # model doit etre un objet résultant de coxph
  tab=matrix(nrow=nrow(summary(model)$coefficients),ncol=ncol(summary(model)$coefficients)+3,
           dimnames = list(c(rownames(summary(model)$coefficients)),c(colnames(summary(model)$coefficients),"RR","IC.inf","IC.sup")))
  tab[,1:ncol(summary(model)$coefficients)]=round(summary(model)$coefficients,digits=3)
  tab[,"RR"]=round(exp(summary(model)$coefficients[,"coef"]),digits=3)
  tab[,"IC.inf"]=round(exp(summary(model)$coefficients[,"coef"]-qnorm(1-alpha/2)*summary(model)$coefficients[,"se(coef)"]),digits=3)
  tab[,"IC.sup"]=round(exp(summary(model)$coefficients[,"coef"]+qnorm(1-alpha/2)*summary(model)$coefficients[,"se(coef)"]),digits=3)
  # attention : en cas de test unilatéral, on divise simplement la p-valeur par 2
  # il faut faire attention au sens du coefficient
  if (sided==1){
    tab[,5]=tab[,5]/2
  }
  # supprimer la colonne exp(coef)
  if (nrow(tab)==1){tab=t(as.matrix(tab[,-2]))} else{tab=tab[,-2]}
  rownames(tab)=rownames(summary(model)$coefficients)
  signif=ifelse(tab[,"Pr(>|z|)"]>=0.1,"",ifelse(tab[,"Pr(>|z|)"]>=0.05,".",ifelse(tab[,"Pr(>|z|)"]>=0.01,"*",ifelse(tab[,"Pr(>|z|)"]>=0.001,"**","***"))))
  res=cbind(tab,signif)
  colnames(res)=c(colnames(tab),"")
  return(noquote(res))
}

#N=100
#time=rexp(N)
#eps=sample(0:1,N,T)
#x=rnorm(N)
#y=rbinom(N,1,0.5)
#
#library(survival)
#IC_RR_coxph(coxph(Surv(time,eps)~x+y))
#IC_RR_coxph(coxph(Surv(time,eps)~x))
