
den01 <- function(x) {dunif(x,min=0,max=1)}
den02 <- function(x) {dexp(x)}
den03 <- function(x) {x*exp(-(x^2)/2)*(x>0)}
den04 <- function(x) {0.5*dexp(x,rate=1)+(x!=0)*0.5*dexp(-x,rate=1)}
den05 <- function(x) {dlogis(x)}
den06 <- function(x) {dcauchy(x,location=0,scale=1)}
den07 <- function(x) exp(-x-exp(-x))
den08 <- function(x) {    outvek<-rep(0,length(x))
  outvek[(x>0)&(x<=1)] <- 1/(2*sqrt(abs(x)+1*(x==0)))[(x>0)&(x<=1)]
  outvek
}
den09 <- function(x) {     outvek<-rep(0,length(x))
  outvek[(x>=1)] <- 1/(2*x[(x>=1)]^(3/2))
  outvek
}
den10 <- function(x) {1/(4*(1+abs(x))^(3/2))}
den11 <- function(x) {dnorm(x,mean=0,sd=1)}
den12 <- function(x) {dlnorm(x, meanlog = 0, sdlog = 1)}
den13 <- function(x) {0.5*dunif(x,min=-0.5,max=0.5) + 0.5*dunif(x,min=-5,max=5)}
den14 <- function(x) { outvek <- 1/(abs(x)*(log(abs(x+1*(x==0)+1*(abs(x)==1)))^2))*(x > -1/(exp(1)^2))*(x < 1/(exp(1)^2))
  outvek[x==0] <- 0
  outvek  
}
den15 <- function(x) {     outvek<-rep(0,length(x))
  outvek[(x>0)&(x<1)] <- -log(x[(x>0)&(x<1)])
  outvek   
}
den16 <- function(x) {(1-abs(x))*((1-abs(x))>0)}
den17 <- function(x) {dbeta(x, shape1=2, shape2=2)}
den18 <- function(x) {dchisq(x,df=1)}
den19 <- function(x) { outvek <- dnorm(sign(x)*abs(x)^(1/3))*(1/3)*((x^2)^(-1/3))
  outvek[x==0] <- 0
  outvek        
}
den20 <- function(x)  {     outvek<-rep(0,length(x))
  outvek[(x>0)] <- exp(-x[(x>0)]^(-1/2))*(0.5*x[(x>0)]^(-3/2))
  outvek
}
den21 <- function(x) {(1/3)*dnorm(x,mean=-20,sd=1/4)+(2/3)*dnorm(x,mean=0,sd=1)}
den22 <- function(x) {(3/4)*dnorm(x,mean=0,sd=1)+(1/4)*dnorm(x,mean=1.5,sd=1/3)}
den23 <- function(x) {(1/2)*dnorm(x,mean=0,sd=1) + (1/10)*dnorm(x,mean=-1,sd=0.1) +(1/10)*dnorm(x,mean=-0.5,sd=0.1)+ 
  (1/10)*dnorm(x,mean=0,sd=0.1) + (1/10)*dnorm(x,mean=0.5,sd=0.1)+(1/10)*dnorm(x,mean=1,sd=0.1)
}
den24 <- function(x) {(32/63)*dnorm(x,mean=-31/21,sd=32/63)+(16/63)*dnorm(x,mean=17/21,sd=16/63)+(8/63)*dnorm(x,mean=41/21,sd=8/63)+
  (4/63)*dnorm(x,mean=53/21,sd=4/63)+(2/63)*dnorm(x,mean=59/21,sd=2/63)+(1/63)*dnorm(x,mean=62/21,sd=1/63)
}
den25 <- function(x) { x2 <- abs(x)+0.1*(abs(x)<0.1)
  outvek <- 2*(1-(x2-0.1)^(1/3))*(abs(x)>=0.1)*(abs(x)<=1.1) 
  outvek
}
den26 <- function(x) {0.5*dunif(x,min=-1,max=1)+0.25*dunif(x,min=20,max=20.1)+0.25*dunif(x,min=-20.1,max=-20)}
den27 <- function(x) { (den16(x+9)+den16(x+7)+den16(x+5)+den16(x+3)+den16(x+1)+den16(x-9)+den16(x-7)+den16(x-5)+den16(x-3)+den16(x-1))/10
}
den28 <- function(x) {
  outvek<-rep(0,length(x))
  outvek[(0<x)&(x<1)]<-(-1/2)*log(abs(x*(1-x))+1*(x==0)+1*(x==1))[(0<x)&(x<1)] 
  outvek
}



`dberdev` <-
function(x,dnum = 1) {
  if (is.nan(dnum) || ! dnum %in% 1:28)
    stop("dnum must be between 1 and 28")
  return (
            eval(               
                parse(text = sprintf("den%02d(x)", dnum)) # evaluate "den[dnum](x)"-string
            )
        )
}
