integrator<-function(init,l,m,psi,times,rtol,atol,migr){
	ode<-function(times,y,p){
		lambda11<-p[1]
		lambda12<-p[2]
		lambda21<-p[3]
		lambda22<-p[4]
		mu1<-p[5]
		mu2<-p[6]
		psi1<-p[7]
		psi2<-p[8]
		
		if (migr==0){
		yd1<-mu1-(lambda11+lambda12+mu1+psi1)*y[1]+lambda11*y[1]*y[1]+lambda12*y[1]*y[2]
		yd2<-mu2-(lambda21+lambda22+mu2+psi2)*y[2]+lambda21*y[1]*y[2]+lambda22*y[2]*y[2]
		yd3<- -(lambda11+lambda12+mu1+psi1)*y[3] + 2*lambda11*y[1]*y[3] + lambda12*y[1]*y[4] + lambda12*y[2]*y[3]
		yd4<- -(lambda22+lambda21+mu2+psi2)*y[4] + 2*lambda22*y[2]*y[4] + lambda21*y[2]*y[3] + lambda21*y[1]*y[4]
		}
		if (migr==1){
		yd1<-mu1-(lambda11+lambda12+mu1+psi1)*y[1]+lambda11*y[1]*y[1]+lambda12*y[2]
		yd2<-mu2-(lambda21+lambda22+mu2+psi2)*y[2]+lambda21*y[1]+lambda22*y[2]*y[2]
		yd3<- -(lambda11+lambda12+mu1+psi1)*y[3] + 2*lambda11*y[1]*y[3] + lambda12*y[4]
		yd4<- -(lambda22+lambda21+mu2+psi2)*y[4] + 2*lambda22*y[2]*y[4] + lambda21*y[3]
		}
		
		#3.4.14
		if (migr==2) {  #SEIR 1 is E, 2 is I
		yd1<- mu1-(lambda12)*y[1]+lambda12*y[2]
		yd2<- mu2-(lambda21+mu2+psi2)*y[2]+lambda21*y[1]*y[2]
		yd3<- -(lambda12)*y[3] + lambda12*y[4]
		yd4<- -(lambda21+mu2+psi2)*y[4] + lambda21*y[1]*y[4]+lambda21*y[2]*y[3]
		}

		
		list(c(yd1,yd2,yd3,yd4))
	}
	out<-lsoda(init,times,ode,c(l,m,psi),rtol=rtol,atol=atol)[2,2:5]
	out
}
