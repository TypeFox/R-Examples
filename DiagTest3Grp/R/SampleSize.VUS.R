SampleSize.VUS <-
function(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,p=0,q=0,lam.minus=1/3,lam0=1/3,lam.plus=1/3,typeIerror=0.05,margin=0.05,subdivisions=50000,...)
  {
    #####This function calculates the sample size to collect data on one biomarker in order to estimate the VUS of the biomarker within margin of error
    ##########Input:
    ###(1)mu.minus,mu0,mu.plus: normal means, can be specified by users arbitrarily (must be in increasing order) or estimated from data first
    ###(2)s.minus,s0,s.plus:normal standard deviations,can be specified by users or estimated from data first
    ###(3)t.minus,t.plus: optimal cut points
    ####(4)p,q: p=minimum desired specificity,and q=minimum desired sensitivity, 0<=p,q<1 ,p,q=0 calculate full VUS, otherwise calculate partial VUS
    ####(5)lam.minus, lam0,lam.plus: for sample size calculation, the expected proportion of samples in the D-, D0 and D+ group, which can be equal or not
    ####(6)typeIerror:type I error rate for sample size calculation,default= 0.05, give 95% CI
    ####(7)margin: for sample size calculation, margin of error on the VUS estimates,,default=0.05. The normal (1-typeIerror)% CI is (VUS-Z_typeIerror*SE(VUS),VUS-Z_typeIerror*SE(VUS)), the sample size calculation will be calculated such that Z_typeIerror*SE(VUS)=margin
    ####(8)subdivisions: # of subintervals for integration using adaptive quadrature
    ####(11)....,  other arguments used in the R function integrate() can be passed along, such as, abs.tol,rel.tol,stop.on.error etc
    
  #Example: SampleSize.VUS(0.5,1.5,3,0.01,0.1,0.5)
    
    #####functions used for integration to obtain VUS estimates, see Xiong et al 2007 paper
    
    f0 <- function(s,a,b,c,d,p,q)
      {
        ###integrate for s over (-Inf, Inf) given a,b,c,d to obtain VUS
        (pnorm(a*s-b)*pnorm(-c*s+d)-p*pnorm(-c*s+d)-q*pnorm(a*s-b)+p*q)*dnorm(s)     
      }

    f1 <- function(s,a,b,c,d,q)
      {
        ##to obtain derivation of VUS w.r.t a
        (dnorm(a*s-b)*pnorm(-c*s+d)-q*dnorm(a*s-b))*dnorm(s)*s
      }


    f2 <- function(s,a,b,c,d,q)
      {
        ##to obtain derivation of VUS w.r.t b
        -(dnorm(a*s-b)*pnorm(-c*s+d)+q*dnorm(a*s-b))*dnorm(s)
      }

    f3 <- function(s,a,b,c,d,p)
      {
        ##to obtain derivation of VUS w.r.t c
        (-pnorm(a*s-b)*dnorm(-c*s+d)+p*dnorm(-c*s+d))*dnorm(s)*s
      }

    f4 <- function(s,a,b,c,d,p)
      {
        ##to obtain derivation of VUS w.r.t c
        (pnorm(a*s-b)*dnorm(-c*s+d)-p*dnorm(-c*s+d))*dnorm(s)
      }

   ###reparametrize the normal sample means and SDs by a,b,c,d
    a <- s0/s.minus
    b <- (mu.minus-mu0)/s.minus
    c <- s0/s.plus
    d <- (mu.plus-mu0)/s.plus


    lower <- (qnorm(p)+b)/a#lower=-Inf when p=q=0
    upper <- (d-qnorm(q))/c
    
    
    V.a <- integrate(f1,a=a,b=b,c=c,d=d,q=q,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
    V.b <- integrate(f2,a=a,b=b,c=c,d=d,q=q,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
    V.c <- integrate(f3,a=a,b=b,c=c,d=d,p=p,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
    V.d <- integrate(f4,a=a,b=b,c=c,d=d,p=p,lower=lower,upper=upper,subdivisions=subdivisions,...)$value
    
    Mv <- 0.5*V.a^2*a^2*(1+lam0/lam.minus)+V.b^2*(a^2+0.5*b^2*lam0/lam.minus+lam0/lam.minus)+0.5*V.c^2*c^2*(1+lam0/lam.plus)+V.d^2*(c^2+0.5*d^2*lam0/lam.plus+lam0/lam.plus)+V.a*V.b*a*b*lam0/lam.minus+V.a*V.c*a*c+2*V.b*V.d*a*c+V.c*V.d*c*d*lam0/lam.plus

    z0 <- qnorm(typeIerror/2,lower.tail=F)
    
    sampleSize <- z0^2*Mv/margin^2

    #sampleSize <- ceiling(sampleSize)
    
    return(sampleSize)
  }

