`power.density.equivalence.md` <-
function(power_sigma,alpha=alpha,theta1=theta1,theta2=theta2,diff=diff,sigma=sigma,n=n,nu=nu){

# power_sigma: x-value for integration   
# alpha:       alpha level for the 2 t-tests (usually 0.05).  alpha level for full test is twice alpha
# theta1:      lower limit of equivalence interval
# theta2:      upper limit of equivalence interval
# diff:        true difference in treatment means for power calculation
# sigma:       root MSE from ANOVA  
# n:           number of observations
# nu:          degrees of freedom for sigma

# integrand for determining power of TOST

pdim<-length(power_sigma)
power_density<-c(1:pdim)
power_t <- qt(1.0-alpha,nu,lower.tail=TRUE, log.p = FALSE)                                                                                                             
  
a<-  sigma/sqrt(nu)    #  constant for transforming chi distribution                                                                                                        
power_const <-  -(nu/2 -1)*log(2) - log(gamma(nu/2)) - log(a)    #  terms not containing power_sigma
                                                                                  
for (isigma in 1:pdim)
{  #  begin for (iterate for the values of power_sigma input by integration program)
  
  power_d1 <-   ((power_sigma[isigma]*power_t*sqrt(2.0/n)) + (theta1-diff))/(sigma*sqrt(2.0/n))                                                                 
  power_d2 <-  ((-power_sigma[isigma]*power_t*sqrt(2.0/n)) + (theta2-diff))/(sigma*sqrt(2.0/n))       
                                                 
  power_phi1 <- pnorm(power_d1,mean=0,sd=1,lower.tail=TRUE, log.p = FALSE)   # terms from difference (normal part)
  power_phi2 <- pnorm(power_d2,mean=0,sd=1,lower.tail=TRUE, log.p = FALSE)                                                                                                      
  power_phi <- power_phi2 - power_phi1   

  power_density[isigma] <- power_phi*exp(power_const + (nu-1)*log(power_sigma[isigma]/(a)) - 0.5*(power_sigma[isigma]/(a))^2 )

}  #  end for

return(power_density)}

