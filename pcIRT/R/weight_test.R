#' Test for the scoring weights in the unidimensional polytomous Rasch model
#' 
#' This functions tests the fit of fixed scoring parameters in a unidimensional
#' polytomous Rasch model.
#' 
#' If the unidimensional polytomous Rasch model fits the data, the weight test
#' can be performed to test whether assumed scoring parameters are appropriate.
#' An unconstrained unidimensional polytomous Rasch model is calculated
#' including estimation of scoring parameters. Furthermore a constrained
#' unidimensional polytomous Rasch model is estimated with fixed scoring
#' parameters (according to the input). Subsequently a Likelihood Ratio test
#' tests the fit of the fixed scoring parameters.
#' 
#' @aliases weight_test summary.wt print.wt
#' @param MPRMobj Object of class \code{MPRM}
#' @param score_param Numerical vector with the scoring parameters that are
#' tested

#' @return \item{emp_Chi2}{\eqn{\chi^2} distributed value of the Likelihood
#' Ratio test} \item{df}{degrees of freedom of the test statistic}
#' \item{pval}{p value of the test statistic}
#' \item{unconstrLoglikelihood}{log-likelihood of the unconstrained model}
#' \item{constrLoglikelihood}{log-likelihood of the constrained model}
#' \item{unconstrNrPar}{number of estimated parameters in the unconstrained
#' model} \item{constrNrPar}{number of estimated parameters in the constrained
#' model} \item{unconstrItempar}{estimated item parameters of the unconstrained
#' model} \item{constrItempar}{estimated item parameters of the constrained
#' model} \item{unconstrScoreParameter}{estimated scoring parameters of the
#' unconstrained model}
#' @author Christine Hohensinn
#' @seealso \code{\link{MPRM}} \code{\link{dLRT}}
#' @references Fischer, G. H. (1974). Einfuehrung in die Theorie
#' psychologischer Tests [Introduction to test theory]. Bern: Huber.
#' @keywords weight test scoring
#' 
#' @rdname wt
#' 
#' @examples
#' 
#' #simulate data set
#' simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), 
#'                   ncol=4),0), 500)
#' 
#' #estimate MPRM item parameters
#' res_mprm <- MPRM(simdat$datmat)
#' 
#' #tests the scoring parameter 0.5 for the unidimensional polytomous model
#' res_weight <- weight_test(res_mprm,  score_param=c(0.5))
#' summary(res_weight)
#' 
#' 
#' @export weight_test
weight_test <-
function(MPRMobj, score_param){
  
  call <- match.call()
  
  if(length(score_param) != (length(table(MPRMobj$data))-2)){stop("Error: wrong number of score parameters!")}
    
  #LR-Test for unidimensionality of item parameters
  
  ep_res  <- EPRM_red(MPRMobj$data)
  ep_resS <- EPRM_red(MPRMobj$data, score_par=score_param)
  
    chi2 <- -2*(ep_resS$logLikelihood-ep_res$logLikelihood)
    df <- length(ep_res$estpar)-length(ep_resS$estpar)
    pvalue <- 1-pchisq(chi2,df)
    
  res <- list(emp_Chi2=chi2, df=df, pval=pvalue, unconstrLoglikelihood=ep_res$logLikelihood, constrLogLikelihood=ep_resS$logLikelihood, unconstrNrPar=length(ep_res$estpar),constrNrPar=length(ep_resS$estpar), unconstrItempar=ep_res$itempar*(-1), constrItempar=ep_resS$itempar*(-1), unconstrScoreParameter=ep_res$score_par)
  class(res) <- "wt"
  res  
}
