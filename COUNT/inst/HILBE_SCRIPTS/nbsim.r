# nbsim.r   
# Table 9.5: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
# Synthetic NB2 Monte Carlo estimation
#
library(MASS)
mysim <- function()
{
 nobs <- 50000
 x1 <-runif(nobs)
 x2 <-runif(nobs)
 xb <- 2 + .75*x1 - 1.25*x2
 a <- .5                        
ia <- 1/.5                     
exb <- exp(xb)                 
xg <- rgamma(n = nobs, shape = a, rate = a)
xbg <-exb*xg                   
nby <- rpois(nobs, xbg)       
nbsim <-glm.nb(nby ~ x1 + x2)  
   alpha <- nbsim$theta 
   pr <- sum(residuals(nbsim, type="pearson")^2)
   prdisp <- pr/nbsim$df.residual
   beta <- nbsim$coef
   list(alpha,prdisp,beta)
}
B <- replicate(100, mysim())   
mean(unlist(B[1,]))
mean(unlist(B[2,]))
apply(matrix(unlist(B[3,]),3,100),1,mean)







