triangular<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
		k<-1-z
		k
		}

epan<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
				z<-q/ban			
				k<-1-z^2
				}
				
bisq<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
		k<-(1-z^2)^2
		}


parzen<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
		k<-z
		tmp1<-which(z<=0.5)
		tmp2<-which(z>0.5)
k[tmp1]<-1-6*z[tmp1]^2+6*abs(z[tmp1])^3
k[tmp2]<-2*(1-abs(z[tmp2]))^3
		k
		}

th<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
k<-(1+cos(pi*z))/2
		k
		}


qs<-function(q,bandwidth){
			if (is.numeric(bandwidth)) ban<-bandwidth
		else ban<-max(q)
		z<-q/ban
k<-(25/(12*pi^2*z^2)) *((sin((6*pi*z)/5)/((6*pi*z)/5))-cos((6*pi*z)/5))
		k
		}
