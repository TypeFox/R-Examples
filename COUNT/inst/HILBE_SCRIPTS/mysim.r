# mysim.r  
# Table 6.5 : Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
#  Monte Carlo simulation - Poisson
mysim <- function()
{
 nobs <- 50000
 x1 <-runif(nobs)
 x2 <-runif(nobs)
py <- rpois(nobs, exp(2 + .75*x1 - 1.25*x2))       
poi <- glm(py ~ x1 + x2, family=poisson)  
   pr <- sum(residuals(poi, type="pearson")^2)
   prdisp <- pr/poi$df.residual
   beta <- poi$coef
   list(beta,prdisp)
}
B <- replicate(100, mysim())   
apply(matrix(unlist(B[1,]),3,100),1,mean)
mean(unlist(B[2,]))
