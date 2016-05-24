qen01 <- function(p) {qunif(p)}
qen02 <- function(p) {qexp(p,rate=1)}    
qen03 <- function(p) {sqrt(-2*log(1-p))}
qen04 <- function(p)  {
  outvek <- rep(0,length(p))
  outvek[p>0.5]<-qexp((abs(2*p-1)),rate=1)[p>0.5]
  outvek[p<0.5]<- -qexp(abs(1-2*p),rate=1)[p<0.5]
  outvek
}       
qen05 <- function(p) {qlogis(p,location=0,scale=1)}
qen06 <- function(p) {qcauchy(p, location = 0, scale = 1)}
qen07 <- function(p) {-log(-log(p))}
qen08 <- function(p) {p^2}
qen09 <- function(p) {1/(1-p)^2}
qen10 <- function(p)  {
  outvek <- rep(0,length(p))        
  outvek[p>0.5]<- (1/(1-abs(2*p-1))^2-1)[p>0.5]
  outvek[p<0.5]<- 1-(1/(1-abs(1-2*p))^2)[p<0.5]
  outvek
}       
qen11 <- function(p) {qnorm(p,mean=0,sd=1)} 
qen12<- function(p) {qlnorm(p,meanlog = 0, sdlog = 1)}
qen13 <- function(p) { (p>=0)*(p<=0.225)*(p*20-5)+(p>0.225)*(p<=0.775)*((p-0.225)/0.55-0.5)+
	(p>0.775)*(p<=1)*((p-0.775)*20+0.5)
}
qen14 <- function(p)  { outvek <- rep(0,length(p))
  outvek[p!=0.5]<-(sign(p-0.5)*(exp(sign(p-0.5)/(0.5-p))))[p!=0.5]
  outvek
}
qen15 <- function(p) { outvek<-0*(p==0)+1*(p==1)
  p2<-unique(p[(p>0)&(p<1)])            
  for (h in p2) {
    fkt<-function(z) (pberdev(z,dnum=15)-h) 
    outvek[p==h]<-uniroot(f=fkt,interval=c(0,1))$root 
  }
  outvek
}
qen16 <- function(p) {(p>=0)*(p<=0.5)*(sqrt(2*p)-1)+(p>0.5)*(p<=1)*(1-sqrt(2-2*p))}
qen17 <- function(p) {qbeta(p, shape1=2, shape2=2)}
qen18 <- function(p) {qchisq(p,df=1)}
qen19 <- function(p) {(qnorm(p,mean=0,sd=1))^3}
qen20 <- function(p) {(p>0)*(1/(log( p*(p>0)+10*(p<0)))^2 )}
qen21 <- function(p) { outvek<-rep(0,length(p))
	for (i in 1:length(p)) {
	p2<-p[i]
	if ((p2>=10^(-5))&(p2<=(1-10^(-5)))) {
	    	fkt<-function(z) (pberdev(z,dnum=21)-p2) 
	    	int<-c(-22,5) 
		}	
	if ((p2>0)&(p2<10^(-5))) {
	    	#left<--21.5
		int<-c(-22.5,-20.5)
            	while (pberdev(int[1],dnum=21)>p2) int<-int-c(1,1)
    	    	fkt<-function(z) (pberdev(z,dnum=21)-p2) 
            	#int<-c(left,-20.5) 
		}
	if ((p2<1)*(p2>(1-10^(-5)))) {
    		#right<-5
		int<-c(4,6)
    		while (pberdev(int[2],dnum=21)<p2) {
      		int<-int+c(1,1)
      		}
	        fkt<-function(z) (pberdev(z,dnum=21)-p2) 
      		}
	if (p2==0) outvek[i]<- -Inf
	if (p2==1) outvek[i]<- Inf
	if ((p2>0)&(p2<1)) outvek[i]<-uniroot(f=fkt,interval=int)$root
	}
  outvek
  }

qen22 <- function(p) { outvek<-rep(0,length(p))
	for (i in 1:length(p)) {
	p2<-p[i]
	if ((p2>=10^(-5))&(p2<=(1-10^(-5)))) {
	    	fkt<-function(z) (pberdev(z,dnum=22)-p2) 
	    	int<-c(-5,5) 
		}	
	if ((p2>0)&(p2<10^(-5))) {
	    	#left<--21.5
		int<-c(-5,-3)
            	while (pberdev(int[1],dnum=22)>p2) int<-int-c(1,1)
    	    	fkt<-function(z) (pberdev(z,dnum=22)-p2) 
            	#int<-c(left,-20.5) 
		}
	if ((p2<1)*(p2>(1-10^(-5)))) {
    		#right<-5
		int<-c(2,4)
    		while (pberdev(int[2],dnum=22)<p2) {
      		int<-int+c(1,1)
      		}
	        fkt<-function(z) (pberdev(z,dnum=22)-p2) 
      		}
	if (p2==0) outvek[i]<- -Inf
	if (p2==1) outvek[i]<- Inf
	if ((p2>0)&(p2<1)) outvek[i]<-uniroot(f=fkt,interval=int)$root
	}
  outvek
  }

qen23 <- function(p) { outvek<-rep(0,length(p))
	for (i in 1:length(p)) {
	p2<-p[i]
	if ((p2>=10^(-5))&(p2<=(1-10^(-5)))) {
	    	fkt<-function(z) (pberdev(z,dnum=23)-p2) 
	    	int<-c(-5,5) 
		}	
	if ((p2>0)&(p2<10^(-5))) {
	    	#left<--21.5
		int<-c(-5,-3)
            	while (pberdev(int[1],dnum=23)>p2) int<-int-c(1,1)
    	    	fkt<-function(z) (pberdev(z,dnum=23)-p2) 
            	#int<-c(left,-20.5) 
		}
	if ((p2<1)*(p2>(1-10^(-5)))) {
    		#right<-5
		int<-c(3,5)
    		while (pberdev(int[2],dnum=23)<p2) {
      		int<-int+c(1,1)
      		}
	        fkt<-function(z) (pberdev(z,dnum=23)-p2) 
      		}
	if (p2==0) outvek[i]<- -Inf
	if (p2==1) outvek[i]<- Inf
	if ((p2>0)&(p2<1)) outvek[i]<-uniroot(f=fkt,interval=int)$root
	}
  outvek
  }


qen24 <- function(p) { outvek<-rep(0,length(p))
	for (i in 1:length(p)) {
	p2<-p[i]
	if ((p2>=10^(-5))&(p2<=(1-10^(-5)))) {
	    	fkt<-function(z) (pberdev(z,dnum=24)-p2) 
	    	int<-c(-5,5) 
		}	
	if ((p2>0)&(p2<10^(-5))) {
	    	#left<--21.5
		int<-c(-5,-3)
            	while (pberdev(int[1],dnum=24)>p2) int<-int-c(1,1)
    	    	fkt<-function(z) (pberdev(z,dnum=24)-p2) 
            	#int<-c(left,-20.5) 
		}
	if ((p2<1)*(p2>(1-10^(-5)))) {
    		#right<-5
		int<-c(3,5)
    		while (pberdev(int[2],dnum=24)<p2) {
      		int<-int+c(1,1)
      		}
	        fkt<-function(z) (pberdev(z,dnum=24)-p2) 
      		}
	if (p2==0) outvek[i]<- -Inf
	if (p2==1) outvek[i]<- Inf
	if ((p2>0)&(p2<1)) outvek[i]<-uniroot(f=fkt,interval=int)$root
	}
  outvek
  }


qen25 <- function(p) {
  outvek <- 1.1*(p==1)-1.1*(p==0)+(-0.1)*(p==0.5)
  p2<-unique(p[(p>0)&(p<1)&(p!=0.5)])           
  for (h in p2) {
    bnd<-c(-1.1,-0.1)*(h<0.5)+c(0.1,1.1)*(h>0.5)
        fkt<-function(z) (pberdev(z,dnum=25)-h) 
        outvek[p==h]<-uniroot(f=fkt,interval=bnd)$root 
        }
  outvek
}
qen26 <- function(p) {(p>=0)*(p<0.25)*(p/2.5-20.1)+(p>=0.25)*(p<0.75)*(4*(p-0.25)-1)+(p>=0.75)*((p-0.75)/2.5+20)} 
qen27<- function(p) {(p<1)*(qen16((10*p)%%1)+c(seq(-9,9,2),10)[floor(10*p)+1]) + 10*(p==1) }
qen28 <- function(p) { outvek<-(-1)*(p==0)+1*(p==1)
  p2<-unique(p[(p>0)&(p<1)])            
  for (h in p2) {
    fkt<-function(z) (pberdev(z,dnum=28)-h) 
        outvek[p==h]<-uniroot(f=fkt,interval=c(-1,1))$root 
        }
  outvek
}

`qberdev` <- function(p, dnum=1) {
        if (is.nan(dnum) || ! dnum %in% 1:28)
            stop("dnum must be between 1 and 28")
        p[p>1]<-1
        p[p<0]<-0
        return (
            eval(               
                parse(text = sprintf("qen%02d(p)", dnum)) # evaluate "qen[dnum](p)"-string
            )
        )
}
