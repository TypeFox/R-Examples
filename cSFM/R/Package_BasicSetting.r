# derivatives and pdf functions for Standard Skewed Normal distribution
# will be merged with "code in the paper.txt"


####################### derivatives used in the model estimation ############
## goal: obtain the 1st and 2nd derivatives of the log-pdf of Y
## Y is SSN (with 0 mean, 1 var, and skewness \gamma)
## i.e. derivatives of log(g(y,gamma)) derivatives wrt y and gamma(skewness) 



### derivatives of standardized SN with mean 0, variance 1 and shape parameter a
## input:   argument y
##          shape parameter a                            
## output:  a list with 6 components (function values and 5 derivatives)  

D.SN <- function(y,a){
  xi<-sqrt(2/pi)*a/sqrt(1+a^2)
  omega<-sqrt(1-xi^2)
  x = omega*y+xi
  
  # pre-requisite function values
  # derivatives using SN when location = 0, scale = 1, shape = a
  pdf.simple <-  2*dnorm(x)*pnorm(x*a)     
  D1.simple <- 1/pi*(-x*pi*pdf.simple + a*sqrt(2*pi)*dnorm(sqrt(1+a^2)*x))
  D11.simple <- -(1 - x^2) * pdf.simple - a * x/pi * (2+a^2) * sqrt(2*pi) * dnorm(sqrt(1+a^2)*x)  
  D2.simple <- x/pi*sqrt(2*pi)*dnorm(sqrt(1+a^2)*x)
  D22.simple <- -1*a*x^2*D2.simple
  D12.simple <- 1/pi*sqrt(2*pi)*dnorm(sqrt(1+a^2)*x)*(1-(1+a^2)*x^2)
  
  # derivatives of SSN with mean 0, variance 1 and shape a
  sn  <- omega*pdf.simple   # function value 
  sn1 <- omega*D1.simple*omega  # 1st derivative wrt. y
  sn11 <- omega^3*D11.simple    # 2nd derivative wrt. y
  
  xi_a<- sqrt(2/pi)*(1+a^2)^(-1.5)
  w_a<- -xi/omega*xi_a
  x_a<- y*w_a + xi_a
  
  #1st derivative wrt. a
  sn2 <- w_a*pdf.simple+ omega*(D1.simple*x_a + D2.simple)  
  # mixture 2nd derivative wrt. y and a
  sn12 <- 2*omega*w_a*D1.simple + omega^2*(D11.simple*x_a + D12.simple)  
  xi_aa = sqrt(2/pi)*(-3)*a*(1+a^2)^(-2.5)
  w_aa = -(1-xi^2)^(-1.5)*xi_a*xi_a - xi/sqrt(1-xi^2)*xi_aa
  x_aa = y*w_aa + xi_aa
  # 2nd derivative wrt. a
  sn22 <- w_aa*pdf.simple + 2*w_a*(D1.simple*x_a+D2.simple) + 
    omega*((D11.simple*x_a + D12.simple) *x_a + 
             D1.simple*x_aa + D12.simple*x_a + D22.simple)
             
  return(list(sn = sn, sn1 = sn1, sn11 = sn11, 
              sn12 = sn12, sn2 = sn2, sn22 = sn22))}

#########################################################
### reparameterization between shape and skewness########
######################################################### 

# transform shape to skewness
skewness.cp <- function(alpha) {        
  delta <- alpha/ sqrt( 1 + alpha^2)
  temp <- delta * sqrt(2/pi)
  (4 - pi)/2 * (temp)^3 / (1 - temp^2)^(3/2)}

# transfrom skewness to shape
shape.dp <- function(gamma)  { 
  b <- sqrt(2/pi)   
  A <- sign(gamma) * (abs(2 * gamma/(4 - pi)))^(1/3)
  delta <- A/(b * sqrt(1 + A^2))
  lambda <- delta/sqrt(1 - delta^2)
  return(lambda)}

#derivatives of gamma wrt shape(alpha)  
D.gamma <- function(alpha){                 
  gamma1 <- -3*pi*(pi-4)*alpha^2*sqrt(2)/ sqrt((pi + pi*alpha^2 - 2*alpha^2)^5)
  
  gamma2 <- 3*(3*pi*alpha^2 - 6*alpha^2 - 2*pi) * pi *(pi - 4) * 
    alpha*sqrt(2) /sqrt((pi + pi*alpha^2 - 2*alpha^2)^7)
  
  return(list(D1 = gamma1, D2 = gamma2))}


##  pdf of  sn(mean = 0, var = 1, skewness = gamma) 
g <-function(y, gamma, log = FALSE) {
  alpha <- shape.dp(gamma)
  delta <- alpha/sqrt(1+alpha^2)
  beta <- sqrt(1-2*delta^2/pi)
  a <- beta*y+sqrt(1-beta^2)*sign(delta)
  
  if (log == FALSE){
    return(2*beta*dnorm(a)*pnorm(alpha*a))
  }  
  if (log == TRUE){
    return(log(2) + log(beta) + dnorm(a, log=T) + pnorm(alpha * a, log.p=T))
  }  
}       


##log(g(y, gamma)) derivatives wrt y and gamma(skewness)  
D.lg <- function(y, gamma){
  fv <- g(y, gamma)
  a <- shape.dp(gamma)  
  D.temp <- D.SN(y,a)    
  D1 <- D.temp$sn1
  D11 <- D.temp$sn11  
  temp <- D.gamma(alpha = a)
  temp1 <- temp[[1]]
  temp2 <- temp[[2]]
  D2 <- D.temp$sn2/temp1
  D12 <- D.temp$sn12/temp1
  D22 <- (D.temp$sn22 - temp2 * D2)/(temp1^2)
  return(list(D1 = D1/fv, D2 = D2/fv, 
              D12 = D12/fv - 1/(fv^2)*D1*D2, D11 = D11/fv - 1/(fv^2)*D1^2, 
              D22 = D22/fv - 1/(fv^2)*D2^2))
}











