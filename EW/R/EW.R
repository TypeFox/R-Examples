##Edgeworth Expansion 
##author:H.R.Law
##email:4islands@gmail.com

##This R-script provides a function EW which generate a 
##polynomial function p_j(x) which approximate the error term 
## of order of o(n^(-j/2)) in F_{W_n} minus Phi(x)
##i.e. We can estimate how fast this standardized Wn converges to 
## N(0,1) using Sigma p_j(x)o(n^(-j/2))
##This result is called Edgeworth Expansion in modern statistical asymptotic theory.

##EW(rvlist,miu,sigma,error,ord)	returns a polynomial function p_ord, assuming the r.v.s holds a mena of miu and a variance of sigma. 
##						error is the eps used to calculate the numerical derivation needed.
##EWplot(f)					plot a graph showing the asymptotic manner of p_ord
##Example:
##source('EW.R')
##data=rnorm(99,0,1)
##fv=EW(data,ord=1)
##EWplot(fv)
EW<-function(rvlist,miu=0,sigma=1,e=10^-5,ord=1){
	rvlist=as.numeric(rvlist)
	##For details of Edgeworth Expansion for error terms in the convergence
	##of F_{W_n} minus Phi(x), refer [Jun Shao]Mathematical Statstics, revised ed, Springer:2003 P70-76, Sec1.5.6
	##n<-length(rvlist)
	##Wn<-sum((rvlist-miu)/sigma)*(n^(-1/2))
	##Cumulant mgf of the sample
	kappa<-function(t){
		log(mgf(t,rvlist))
	}
	##Generate the moment generating function based on the sample
	mgf<-function(t,rv=c()){
		mean(exp(t%*%rv))
	}
	##Numerical differentiation
	diff<-function(f,e=10^-5){
		ft<-function(x){
			(f(x+e)-f(x-e))/(2*e)
		}
		ft
	}
	kappa3<-diff(diff(diff(kappa)))
	kappa4<-diff(diff(diff(diff(kappa))))
	##Calculate the kappa value using the numerical derivative at zero
	k3<-kappa3(0)*3*2*1*e^3
	k4<-kappa4(0)*4*3*2*1*e^4
	##This is the derivative of the cdf of N(0,1), I deal with it separatly
	PhiP<-function(tkk,e=10^-5){
		##Define the cdf of N(0,1)
		Phi<-function(val1){
			pnorm(val1,0,1)
		}
		(Phi(tkk+e)-Phi(tkk-e))/(2*e)
	}
	##Calculate the first two polynomials in the Edgeworth Expansion
	p1<-function(y){
		(-1/6)*k3*(y^2-1)*PhiP(y)
	}
	p2<-function(y){
		-((-1/24)*k4*y*(y^2-3)+(1/72)*k3*y*(y^4-10*y^2+15))*PhiP(y)
	}
	if(ord==1)p1
	else p2
	##p1(0)
	##p2(0)
}
##Plotting function for Edgeworth polynomials
EWplot<-function(f){
	v<-seq(-5,5,by=.2)
	plot(v,f(v),font.lab=2,type="o")
}