pen01 <- function(q) {punif(q)}
pen02 <- function(q) {pexp(q,rate=1)}    
pen03 <- function(q) {(q>0)*(1-exp(-q^2/2))}        
pen04 <- function(q) {(q<0)*0.5*(1-pexp(-q,rate=1))+(q>0)*0.5*pexp(q,rate=1)+0.5*(q>=0)}
pen05 <- function(q) {plogis(q,location=0,scale=1)}
pen06 <- function(q) {pcauchy(q, location = 0, scale = 1)}   
pen07 <- function(q) {exp(-exp(-q))}
pen08 <- function(q) {(q>0)*(q<1)*sqrt(abs(q))+1*(q>=1)}
pen09 <- function(q) {(q>1)*(1-1/sqrt(abs(q)+1*(q==0)))}
pen10 <- function(q) {(q>=0)*(1-1/(2*sqrt(abs(q)+1)))+(q<0)*(1/(2*sqrt(1+abs(q))))} 
pen11 <- function(q) {pnorm(q,mean=0,sd=1)} 
pen12 <- function(q) {plnorm(q,meanlog = 0, sdlog = 1)}
pen13 <- function(q) {0.5*punif(q,min=-0.5,max=0.5)+0.5*punif(q,min=-5,max=5)}
pen14 <- function(q) {   outvek <- rep(0,length(q))
  outvek[q==0]<-0.5
  outvek[(q!=0)&(abs(q)<exp(-2))] <- (0.5-sign(q)/log(abs(q)))[(q!=0)&(abs(q)<exp(-2))]
  outvek[q>exp(-2)]<-1
  outvek[q==Inf]<-1
  outvek[q==-Inf]<-0
  outvek
}
pen15 <- function(q) { outvek<-1*(q>=1)
  outvek[(q>0)&(q<1)] <- (q-q*log(abs(q)))[(q>0)&(q<1)]
  outvek[q==Inf]<-1
  outvek[q==-Inf]<-0
  outvek        
}
pen16 <- function(q) {(q>(-1))*(q<0)*(0.5+q+0.5*q^2)+(q>=0)*(q<1)*(0.5+q-0.5*q^2)+(q>=1)}
pen17 <- function(q) {pbeta(q, shape1=2, shape2=2)}
pen18 <- function(q) {pchisq(q,df=1)}
pen19 <- function(q) {pnorm(sign(q)*abs(q)^(1/3),mean=0,sd=1)}
pen20 <- function(q) {(q>0)*exp(-1/sqrt( q*(q>0)+10*(q<=0)  ))}
pen21 <- function(q) {(1/3)*pnorm(q,mean=-20,sd=1/4)+(2/3)*pnorm(q,mean=0,sd=1)}
pen22 <- function(q) {(3/4)*pnorm(q,mean=0,sd=1)+(1/4)*pnorm(q,mean=1.5,sd=1/3)}
pen23 <- function(q) {0.5*pnorm(q,mean=0,sd=1)+0.1*pnorm(q,mean=-1,sd=0.1)+0.1*pnorm(q,mean=-0.5,sd=0.1)+0.1*pnorm(q,mean=0,sd=0.1)+
  0.1*pnorm(q,mean=0.5,sd=0.1)+0.1*pnorm(q,mean=1,sd=0.1)
}
pen24 <- function(q) {(32/63)*pnorm(q,mean=-31/21,sd=32/63)+(16/63)*pnorm(q,mean=17/21,sd=16/63)+(8/63)*pnorm(q,mean=41/21,sd=8/63)+
  (4/63)*pnorm(q,mean=53/21,sd=4/63)+(2/63)*pnorm(q,mean=59/21,sd=2/63)+(1/63)*pnorm(q,mean=62/21,sd=1/63)
}
pen25 <- function(q) {0.5*(q>=1.1)+(q>=-0.1)*0.5+(q<1.1)*(q>0.1)*0.5*(4*(q-0.1)-3*(abs(q-0.1)^(4/3)))+
  (q>-1.1)*(q<(-0.1))*0.5*(1-(4*(abs(q)-0.1)-3*abs(-q-0.1)^(4/3)))
}
pen26 <- function(q) {0.5*punif(q,min=-1,max=1)+0.25*punif(q,min=20,max=20.1)+0.25*punif(q,min=-20.1,max=-20)}
pen27 <- function(q){ (pen16(q+9)+pen16(q+7)+pen16(q+5)+pen16(q+3)+pen16(q+1)+pen16(q-9)+pen16(q-7)+pen16(q-5)+pen16(q-3)+pen16(q-1))/10
}
pen28 <- function(q)    {outvek<-(q>=1)
  outvek[(q>0)&(q<1)]<-(-0.5*q*log(abs(q))+0.5*(1-q)*log(abs(1-q))+q)[(q>0)&(q<1)]
  outvek[q==Inf]<-1
  outvek[q==-Inf]<-0
  outvek
}

`pberdev` <- function(q, dnum=1) {
        if (is.nan(dnum) || ! dnum %in% 1:28)
            stop("dnum must be between 1 and 28")
        return (
            eval(               
                parse(text = sprintf("pen%02d(q)", dnum)) # evaluate "pen[dnum](q)"-string
            )
        )
}
