###################################################################################
### Purpose: Calulate density given hyperparemeters and model
### input: hyperparameters(a1,b1,a2,b2), data(y1,n1,y2,n2), other(measure,model),
###        theta (comparative measure)
### Output: density(at theta)
### Note: the function uses functions  "dens.sar" and "dens.inde"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
densityfunction <- function(theta=theta,a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,y1=0,y2=0,n1=0,
                            n2=0,measure=measure,model=model) {
  if(model=="Sarmanov") {
    cc <- sqrt(a1*a2*b1*b2)/sqrt((a1+b1+1)*(a2+b2+1))
    upper.bound <- cc/max(a1*b2, a2*b1)
    lower.bound <- -cc/max(a1*a2, b1*b2)
    if (rho > upper.bound | rho < lower.bound) stop(paste("rho is out of bound: ",
                                                          lower.bound, upper.bound))
  }
  alpha1 <- y1+a1; beta1 <- n1-y1+b1
  alpha2 <- y2+a2; beta2 <- n2-y2+b2

  if (model=="Sarmanov") result <- dens.sar(theta=theta,a1=a1,b1=b1,a2=a2,b2=b2,
        rho=rho,y1=y1,y2=y2,n1=n1,n2=n2,measure=measure)
  if (model=="Independent") result <- dens.inde(alpha1,beta1,alpha2,beta2,theta,measure)      
  return(result)
}

###################################################################################
### Purpose: Calulate density in an interval given hyperparemeters
### input: hyperparameters(a1,b1,a2,b2), data(y1,n1,y2,n2), measure, model,
###        interval(begin,end,length)
### Output: list, containning x and correspondent density y
### Note: the function uses functions "dens.sar" and "dens.inde"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
dens.post <- function(y1=y1, n1=n2, y2=y2, n2=n2, a1=a1, b1=b1, a2=a2, b2=b2, rho=rho,
                      grid.start=grid.start, grid.end=grid.end, grid.num=grid.num,
                      measure=measure,model=model) {
  theta.grid <- seq(grid.start, grid.end, length=grid.num)
  theta.dens.prior <- rep(NA, length=length(theta.grid))
  
  for(i in 1:length(theta.grid)) {
    theta.dens.prior[i] <- densityfunction(theta=theta.grid[i],a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,
                                           y1=y1,y2=y2,n1=n1,n2=n2,measure=measure,model=model)
  }
  theta.dens.prior <- as.numeric(theta.dens.prior)
  return(list(y=theta.dens.prior, x=theta.grid))
}

###################################################################################
### Purpose: Calulate the hypergeometric function 2f1 using Fortran library
### input: the parameters required by hypergeometric function 2f1
### Output: the value of hypergeoFun at X=x
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
hypergeoFun <- function(aa, bb, cc, xx, YY=0) {

if(.Platform$OS.type=="unix")
	path<-file.path(find.package("mmeta"),"libs","mmeta.so")
	
  if(.Platform$OS.type=="windows"){
	if(.Platform$r_arch=="i386") path<-file.path(find.package("mmeta"),"libs","i386","mmeta.dll")
	if(.Platform$r_arch=="x64")  path<-file.path(find.package("mmeta"),"libs","x64","mmeta.dll")
  
  }
    dyn.load(path)
	.Fortran("hygfx", a=as.double(aa),b=as.double(bb),c=as.double(cc),
           x=as.double(xx),y=as.double(YY),PACKAGE="mmeta")$y
 
}

###################################################################################
### Purpose: compute the omegas used in the posterior distribution of p1 and p2
###          it is identical for OR/RR/RD.
### input: data (y1, n1, y2, n2), hyperparameters(a1,b1,a2,b2,rho)
### Output: the value of hypergeoFun at X=x
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
omega.computation <- function(y1, n1, y2, n2, a1, b1, a2, b2, rho) {
  alpha1 <- y1+a1; beta1 <- n1-y1+b1
  alpha2 <- y2+a2; beta2 <- n2-y2+b2

  # mu1, mu2: marginal means of p1 and p2
  mu1 <- a1/(a1+b1); mu1.1 <- 1-mu1
  mu2 <- a2/(a2+b2); mu2.1 <- 1-mu2

  # delta1, delta2: marginal sd of p1 and p2
  delta1 <- sqrt(mu1*mu1.1/(a1+b1+1))
  delta2 <- sqrt(mu2*mu2.1/(a2+b2+1))

  # myd: d=(mu1*mu2)/(delta1*delta2)
  myd <- (mu1*mu2)/(delta1*delta2)
  
  ## v1-v4 are weights
  v2 <- v3 <- -rho*myd
  v1 <- 1 - v2
  v4 <- -v2

  temp1 <- (lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2)+lgamma(beta2)
            +lgamma(a1+b1)+lgamma(a2+b2)
           -(lgamma(a1)+lgamma(b1)+lgamma(a2)+lgamma(b2)
            +lgamma(alpha1+beta1)+lgamma(alpha2+beta2)))
  temp2 <- (lgamma(alpha1+1)+lgamma(beta1)+lgamma(alpha2)+lgamma(beta2)
            +lgamma(a1+b1+1)+lgamma(a2+b2)
           -(lgamma(a1+1)+lgamma(b1)+lgamma(a2)+lgamma(b2)
            +lgamma(alpha1+beta1+1)+lgamma(alpha2+beta2)))
  temp3 <- (lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2+1)+lgamma(beta2)
            +lgamma(a1+b1)+lgamma(a2+b2+1)
           -(lgamma(a1)+lgamma(b1)+lgamma(a2+1)+lgamma(b2)
            +lgamma(alpha1+beta1)+lgamma(alpha2+beta2+1)))
  temp4 <- (lgamma(alpha1+1)+lgamma(beta1)+lgamma(alpha2+1)+lgamma(beta2)
            +lgamma(a1+b1+1)+lgamma(a2+b2+1)
           -(lgamma(a1+1)+lgamma(b1)+lgamma(a2+1)+lgamma(b2)
            +lgamma(alpha1+beta1+1)+lgamma(alpha2+beta2+1)))

  eps <- 1e-30
  if (abs(v1) < eps) {
    omega1 <- 0
  } else {
    omega1 <- 1/(1 + v2/v1*exp(temp2-temp1) + v3/v1*exp(temp3-temp1) + v4/v1*exp(temp4-temp1))
  }
  if (abs(v2) < eps) {
    omega2 <- 0
  } else {
    omega2 <- 1/(v1/v2*exp(temp1-temp2) + 1 + v3/v2*exp(temp3-temp2) + v4/v2*exp(temp4-temp2))
  }
  if (abs(v3) < eps) {
    omega3 <- 0
  } else {
    omega3 <- 1/(v1/v3*exp(temp1-temp3) + v2/v3*exp(temp2-temp3) + 1 + v4/v3*exp(temp4-temp3))
  }
  if (abs(v4) < eps) {
    omega4 <- 0
  } else {
    omega4 <- 1/(v1/v4*exp(temp1-temp4) + v2/v4*exp(temp2-temp4) + v3/v4*exp(temp3-temp4) + 1)
  }
  
  return(list(omega1=omega1, omega2=omega2, omega3=omega3, omega4=omega4))
}

###################################################################################
### Purpose: compute the densty of OR for independent model
### input: alpha1, beta1, alpha2, beta2, theta
### Output: the density at theta
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
OR.dens.inde <- function(alpha1, beta1, alpha2, beta2, theta) {
  if(theta>=1) {
    mylog <- ((-1-beta2)*log(theta)
              -(lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2)+lgamma(beta2)+lgamma(alpha1+beta1+alpha2+beta2))
              +(lgamma(alpha1+beta1)+lgamma(alpha2+beta2)+lgamma(alpha1+alpha2)+lgamma(beta1+beta2))
              +log(hypergeoFun(alpha2+beta2, beta1+beta2, alpha1+alpha2+beta1+beta2, 1-(1/theta))))
    result <- exp(mylog)
  }
  if(theta<1) {
    mylog <- ((-1+alpha2)*log(theta)
              -(lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2)+lgamma(beta2)+lgamma(alpha1+beta1+alpha2+beta2))
              +(lgamma(alpha1+beta1)+lgamma(alpha2+beta2)+lgamma(alpha1+alpha2)+lgamma(beta1+beta2))
              +log(hypergeoFun(alpha2+beta2, alpha1+alpha2, alpha1+alpha2+beta1+beta2, 1-theta)))    
    result <- exp(mylog)
  }
  return(result)
}

###################################################################################
### Purpose: compute the densty of RR for independent model
### input: alpha1, beta1, alpha2, beta2, theta
### Output: the density at theta
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
RR.dens.inde <- function(alpha1, beta1, alpha2, beta2, theta){            
  if(theta < 1){
    mylog <- ((-1+alpha2)*log(theta)
              -(lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2)
                +lgamma(beta2)+lgamma(alpha1+alpha2+beta1))
              +(lgamma(alpha1+beta1)+lgamma(alpha2+beta2)
                +lgamma(alpha1+alpha2)+lgamma(beta1))
              +log(hypergeoFun(1-beta2,alpha1+alpha2,alpha1+alpha2+beta1,theta)))
    result <- exp(mylog)
  }
  if(theta >= 1){
    mylog <- ((-1-alpha1)*log(theta)
              -(lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2)
                +lgamma(beta2)+lgamma(alpha1+alpha2+beta2))
              +(lgamma(alpha1+beta1)+lgamma(alpha2+beta2)
                +lgamma(alpha1+alpha2)+lgamma(beta2))
              +log(hypergeoFun(1-beta1,alpha1+alpha2,alpha1+alpha2+beta2,1/theta)))
    result <- exp(mylog)
  }
  return(result)
}


###################################################################################
### Purpose: compute the densty of RR/RR for independent model
### input: alpha1, beta1, alpha2, beta2, theta, meansure
### Output: the density at theta
### Note: this function is used by "OR.dens.inde","RR.dens.inde"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
dens.inde <- function(alpha1,beta1,alpha2,beta2,theta,measure) {
 if (measure=="OR")  result <- OR.dens.inde(alpha1,beta1,alpha2,beta2,theta)
 if (measure=="RR")  result <- RR.dens.inde(alpha1,beta1,alpha2,beta2,theta)

 return(result)
}

###################################################################################
### Purpose: compute the densty of OR under Sarmanov beta distribution
### input: alpha1, beta1, alpha2, beta2, theta, omega1, omega2, omega3, omega4
### Output: the density at theta
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
OR.dens.sar <- function(alpha1,beta1,alpha2,beta2,theta,
                        omega1,omega2,omega3,omega4) {
  return(omega1*OR.dens.inde(alpha1,beta1,alpha2,beta2,theta)
         +omega2*OR.dens.inde(alpha1+1,beta1,alpha2,beta2,theta)
         +omega3*OR.dens.inde(alpha1,beta1,alpha2+1,beta2,theta)
         +omega4*OR.dens.inde(alpha1+1,beta1,alpha2+1,beta2,theta))                 
}

###################################################################################
### Purpose: compute the densty of RR under Sarmanov beta distribution
### input: alpha1, beta1, alpha2, beta2, theta, omega1, omega2, omega3, omega4
### Output: the density at theta
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
RR.dens.sar <- function(alpha1,beta1,alpha2,beta2,theta,
                        omega1,omega2,omega3,omega4) {
  results <- (omega1*RR.dens.inde (alpha1,beta1,alpha2,beta2,theta)
              +omega2*RR.dens.inde (alpha1+1,beta1,alpha2,beta2,theta)            
              +omega3*RR.dens.inde( alpha1,beta1,alpha2+1,beta2,theta)            
              +omega4*RR.dens.inde (alpha1+1,beta1,alpha2+1,beta2,theta))
  return(results)
}

#########################################################################################
### Purpose: compute the densty of OR and RR under Sarmanov beta distribution
### input: hyperparameters (a1, b1, a2, b2, rho), data (y1, n1, y2, n2), measure
###        theta
### Output: the density at theta
### Note: The function uses functions "OR.dens.sar","RR.dens.sar"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
#########################################################################################
dens.sar <- function(theta=theta,a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,y1=y1,y2=y2,n1=n1,n2=n2,
                     measure=measure) {                                                     
  alpha1 <- y1+a1; beta1 <- n1-y1+b1
  alpha2 <- y2+a2; beta2 <- n2-y2+b2
  myOmega <- omega.computation(y1, n1, y2, n2, a1, b1, a2, b2, rho)
  omega1 <- myOmega$omega1
  omega2 <- myOmega$omega2
  omega3 <- myOmega$omega3
  omega4 <- myOmega$omega4
  if (measure=="OR")  {result <- OR.dens.sar(alpha1,beta1,alpha2,beta2,theta,omega1,omega2,omega3,omega4)}
  if (measure=="RR")  {result <- RR.dens.sar(alpha1,beta1,alpha2,beta2,theta,omega1,omega2,omega3,omega4)}
  return(result)
}
