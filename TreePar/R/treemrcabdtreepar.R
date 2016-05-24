treemrcabdtreepar <-
function(x,l,mu,rho) {
	#x<-sort(x)
	x<-c(NA,x)
	a=mu/l
	r=l-mu
	N<-length(x)
    res <- -((N - 2) * log(r*rho) + 
            N * log(1 - a) + r * sum(x[2:(N-1)])  -
             2 * sum(log(rho*exp(r * x[2:N]) +((1-rho)-a))))
    print("this function conditions on survival. shifts and densdep do not!")
    res
}

