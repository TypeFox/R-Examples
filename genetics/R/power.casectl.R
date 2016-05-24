# power calculation for case-control design (default: case:control = 1:1)
# Author: Michael Man
#   Date: May 5, 2004
#      N: total number of subjects
#  gamma: relative risk in multiplicative model;
#         not used in Dominant or Recessive model (assume A as protective allele)
#      p: frequency of A allele
#     kp: prevalence of disease
#  alpha: significance level
#     fc: fraction of cases
#     pi: probability of 'aa' genotype has the disease
#   minh: mode of inheritance
# reference: Long, A. D. and C. H. Langley (1997). Genetic analysis of complex traits. Science 275: 1328.
#            Agresti, A. (2002) Categorical Data Analysis. Second Edition, p243.
# ( modified from pbsize{gap} )
# requirement: It is recommended to use R 1.9.0 or above.  The function 'qchisq' in earlier version
#              has problem with large noncentrality parameter.


# under HWE                       AA       Aa       aa       
#   fHW = p(genotype)        = c( p^2,     2pq,     q^2 )     

# model specification 
#   f.mod = relative risk    = c(gamma^2,  gamma,    1 )    # multiplicative model
#   f.mod =                  = c(  0,        0,      1 )    #       dominant model
#   f.mod =                  = c(  0,        1,      1 )    #      recessive model

# conditional prob.
#   p(D|genotype) = f.mod*pi = c(gamma^2,  gamma,    1 )*pi

# population joint prob. (f.mod = 1 under Ho)
#   Case     p(D,     genotype) = p(genotype)*     p(D|genotype)  = fHW*   f.mod*pi
#   Control  p(D_not, genotype) = p(genotype)*(1 - p(D|genotype)) = fHW*(1-f.mod*pi)

# population conditional prob. (f.mod = 1 under Ho)
#   Case     p(genotype|D)     = p(D    , genotype)/P(D    ) = P(D    , genotype)/sum(P(D    , genotype)) = fHW*   f.mod*pi  /    sum(fHW*f.mod*pi)
#   Control  p(genotype|D_not) = p(D_not, genotype)/P(D_not) = P(D_not, genotype)/sum(P(D_not, genotype)) = fHW*(1-f.mod*pi) / (1-sum(fHW*f.mod*pi))

# sample or allocation probability
#   1:1 case-control design  p(D|Sample) = fc = 1/2
#   1:2 case-control design                fc = 1/3
#   a prospective design                   fc = sum(fHW*f.mod*pi)

# sample joint prob. (f.mod = 1 under Ho)
# for prospective design, this is the same as population joint prob. since 'fc' cancels out with 'sum(fHW*f.mod*pi)' 
#   Case     p(genotype,D    |sample) = p(genotype|D    )*     p(D|Sample)  =    fc *fHW*   f.mod*pi  /    sum(fHW*f.mod*pi)
#   Control  p(genotype,D_not|sample) = p(genotype|D_not)*(1 - p(D|Sample)) = (1-fc)*fHW*(1-f.mod*pi) / (1-sum(fHW*f.mod*pi))


## power.casectrl <- function (N, gamma = 4.5, p = 0.15, kp=.1, alpha=.05, fc=0.5,
##                              minh=c('multiplicative', 'dominant','recessive','partialrecessive')) 
## {
##   minh <- match.arg(minh)
##   if ( !all(gamma > 0, N > 0) ) stop('N and gamma must be greater than 0')
##   if ( min(p, kp, alpha, fc) <= 0 | max(p, kp, alpha, fc) >=1 ) stop('p, kp, alpha, and fc must be between 0 and 1.') 
##   f.mod <- switch(minh,
##          multiplicative = c(gamma^2, gamma, 1),
##       partialrecessive = c(gamma,       1, 1),           
##          dominant       = c(      0,     0, 1),
##          recessive      = c(      0,     1, 1)  ) 
##   q <- 1 - p
##   fhw <- c(p^2, 2*p*q, q^2)
##   pi <- kp/sum(f.mod*fhw)
##   if (pi <= 0 | pi >=1) {
##     warning('The combination of p, kp, and gamma produces an unrealistic value of pi.')
##     ret <- NA
##   } else {
##     fe  <- rbind(fhw, fhw)
##     dimnames(fe) <- list(c("Case", "Control"), c("AA", "Aa", "aa"))
##     f <- fe*rbind(f.mod*pi, 1-f.mod*pi)
##     Pct <- apply(f, 1, sum)
##     f2 <-  f *c(fc, 1-fc)/Pct   # normalize the frequencies for each row
##     fe2 <- fe*c(fc, 1-fc)  
##     fe2; apply(fe2, 1, sum); f2; apply(f2, 1, sum)  
##     lambda <- sum((f2-fe2)^2/fe2)*N
##     ret <- 1 - pchisq(qchisq(1-alpha, df=1), df=1, ncp=lambda, lower.tail=T)
##   }
##   ret
## }

power.casectrl <- function (...)
{
  .Deprecated("'GPC', 'GeneticPower.Quantitative.Factor', or 'GeneticPower.Quantitative.Numeric in the BioConductor GeneticsDesign package")
}
