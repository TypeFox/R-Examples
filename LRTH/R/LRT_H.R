#' The Function for Likelihood Ratio Test Accounting for Genetic Heterogeneity
#' 
#' It gives the asymptotic p-value of the LRT_H test.
#' @param x a n x 1 vector of genotypic score for SNP (i.e. 0, 1 or 2, the number of minor alleles of a SNP); n is the number of observations.
#' @param y a n x 1 vector of disease status; case/xcontrol: 1/0; ; n is the number of observations.
#' @return The asymptotic p-value of LRT_H test.
#' @details Missing values in either x or y (i.e. genotype or disease status) will be removed.
#' @author Zhiyuan (Jason) Xu and Wei Pan
#' @references 
#' Qian M., Shao Y., 2013. A Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity. 
#' Annals of Human Genetics, 77(2): 174-182.
#' 
#' Zhou H., Pan W., 2009. Binomial Mixture Model-based Association Tests under Genetic Heterogeneity. 
#' Annals of Human Genetics, 73(6): 614-630.
#' @examples
#'        y = c(rep(1,500),rep(0,500))
#'        x1 = sample(c(0,1,2),500,replace=TRUE,prob = c(0.64,0.32,0))
#'        x2 = sample(c(0,1,2),500,replace=TRUE,prob = c(0.49,0.42,0))
#'        x = c(x1,x2)
#'        LRT_H(x,y)
#' @export
#' @importFrom stats complete.cases pchisq


LRT_H = function(x,y){
  ## implement Qian & Shao (2013) LRT test
  ## given dat including: x (genotype: 0,1,2) and y(Case/control:1/0)
  ## output p-value
  dat = data.frame(y=y,x=x)
  dat = as.data.frame(dat)
  dat = dat[complete.cases(dat),]
  n0 = nrow(subset(dat,y==1 & x==0))
  n1 = nrow(subset(dat,y==1 & x==1))
  n2 = nrow(subset(dat,y==1 & x==2))
  m0 = nrow(subset(dat,y==0 & x==0))
  m1 = nrow(subset(dat,y==0 & x==1))
  m2 = nrow(subset(dat,y==0 & x==2))
  n = n0+n1+n2
  m = m0+m1+m2
  
  p0 = (n2+m2+(n1+m1)/2)/(n+m)
  
  logL0 = (m0+n0)*log((1-p0)^2) + (m1+n1)*log(2*p0*(1-p0)) + (m2+n2)*log(p0^2)
  ph = (m2+m1/2)/m
  
  logLh = (m0)*log((1-ph)^2) + (m1)*log(2*ph*(1-ph)) + (m2)*log(ph^2)
  
  pd = (n2+n1/2)/n
  
  
  logLd1 = (n0)*log(n0/n) + (n1)*log(n1/n) + (n2)*log(n2/n)
  logLd2 = (n0)*log((1-pd)^2) + (n1)*log(2*pd*(1-pd)) + (n2)*log(pd^2)
  if( (4*n0*n2) > (n1^2) ) logLd = logLd1 else logLd = logLd2
  
  
  twolambda = 2*(logLd+logLh-logL0)
  
  pvalue = (pchisq(twolambda,1,lower.tail=F) + pchisq(twolambda,2,lower.tail=F))/2
  
  pvalue
}

