#' Dimensionality test for the multidimensional polytomous Rasch model
#' 
#' This function tests whether the multidimensional polytomous Rasch model can
#' be reduced to a unidimensional polytomous model.
#' 
#' For this test, a unidimensional model assuming the categories as linearly
#' dependent is computed. Subsequently a Likelihood Ratio test is conducted.
#' 
#' @aliases dLRT summary.dLR print.dLR
#' @param MPRMobj Object of class \code{MPRM}
#' @return \item{emp_Chi2}{\eqn{\chi^2} distributed value of the Likelihood
#' Ratio test} \item{df}{degrees of freedom of the test statistic}
#' \item{pval}{p value of the test statistic}
#' @author Christine Hohensinn
#' @seealso \code{\link{MPRM}} \code{\link{LRT}}
#' @references Fischer, G. H. (1974). Einfuehrung in die Theorie
#' psychologischer Tests [Introduction to test theory]. Bern: Huber.
#' @keywords dimensionality model test
#' 
#' @export
#' @rdname dLR
#' @examples
#' 
#' #simulate data set
#' simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), 
#'    ncol=4),0), 500)
#' 
#' #estimate MPRM item parameters
#' res_mprm <- MPRM(simdat$datmat)
#' 
#' res_dlrt <- dLRT(res_mprm)
#' summary(res_dlrt)
#' 
#' 
#' @export dLRT
dLRT <-
function(MPRMobj){
  
  eprm <- EPRM_red(MPRMobj$data)
  
  emp_Chi2 <- -2*(eprm$logLikelihood-MPRMobj$logLikelihood)
  df <- length(MPRMobj$estpar)-length(eprm$estpar)
  pvalue <- 1-pchisq(emp_Chi2,df)
  
  res <- list(emp_Chi2=emp_Chi2, df=df, pvalue=pvalue)
  class(res) <- "dLR"
  res
}
