cAnsBrad <-
function(alpha,m,n,method=NA,n.mc=10000){

pcalc <- function(q, m, n) .Call("pAnsari", q, m, n)
ccalc <- function(q, m, n) .Call("qAnsari", q, m, n)

N=m+n

if(alpha>1||alpha<0||class(alpha)!="numeric"){
	cat('Error: Check alpha value! \n')
	return(alpha)
}
outp<-list()
outp$m<-m
outp$n<-n
outp$alpha<-alpha
outp$stat.name<-"Ansari-Bradley C"


outp$method<-method
if(is.na(outp$method)){
    if(outp$m+outp$n<=200){
      outp$method<-"Exact"
    }
    if(outp$m+outp$n>200){
      outp$method<-"Asymptotic"
    }
}

outp$two.sided<-1;

if(outp$method=="Monte Carlo"){
	warning("The exact computation will work for large data, so Monte Carlo methods
		are not recommended for this procedure.")
  outp$method="Exact"
}

if(outp$method=="Exact"){
	if(N%%2==0) tot<-N/2*(N/2+1)
	if(N%%2==1) tot<-(N-1)/2*((N-1)/2+1)+(N+1)/2

	outp$cutoff.U<-tot-ccalc(alpha,m,n)+1
	outp$cutoff.L<-tot-ccalc(1-alpha,m,n)-1
	outp$true.alpha.U<-pcalc(tot-outp$cutoff.U,m,n)
	outp$true.alpha.L<-1-pcalc(tot-outp$cutoff.L-1,m,n)
}


if(outp$method=="Asymptotic"){

	if(N%%2==0){
		exp_C=n*(N+2)/4
		var_C=m*n*(N+2)*(N-2)/(48*(N-1))
	}

	if(N%%2==1){
		exp_C=n*(N+1)^2/(4*N)
		var_C=m*n*(N+1)*(3+N^2)/(48*N^2)
	}
	outp$cutoff.L=exp_C+qnorm(alpha)*(var_C)^(1/2)
	outp$cutoff.U=exp_C+qnorm(1-alpha)*(var_C)^(1/2)
}


class(outp)<-"NSM3Ch5c"
outp
}
