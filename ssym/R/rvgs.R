rvgs <-
function(n, family, xi){

if(family!="Normal" & family!="Slash" & family!="Hyperbolic" & family!="Sinh-t" & family!="Sinh-normal" & family!="Contnormal" & family!="Powerexp" & family!="Student")
stop("family of distributions specified by the user is not supported!!",call.=FALSE)

if(family=="Slash" | family=="Hyperbolic" | family=="Sinh-t" | family=="Sinh-normal" | family=="Contnormal" | family=="Powerexp" | family=="Student"){
  if(missingArg(xi)) stop("for the family of distributions specified by the user an extra parameter value is required!!", call.=FALSE) 
}
if(n!=floor(n) | n<=0) stop("n must be a positive integer!!",call.=FALSE)

	if(family=="Normal"){
		rvg <- function(n){
			   rnorm(n)
		}
	}
	if(family=="Student"){
		if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
		rvg <- function(n){
		      rnorm(n)/sqrt(rgamma(n,shape=(xi[1]/2),scale=(2/xi[1])))
		}
	}
	if(family=="Contnormal"){
		if(xi[1]<=0  | xi[1]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
		if(xi[2]<=0  | xi[2]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
		rvg <- function(n){
			   u <- runif(n)
               sqrt(1/ifelse(u<=xi[1],xi[2],1))*rnorm(n)
		}
	}

	if(family=="Powerexp"){
		if(xi[1]<=-1  | xi[1]>=1) stop("the extra parameter must be within the interval (-1, 1)!!",call.=FALSE)
		rvg <- function(n){
			p <- 2/(xi[1]+1)
	        sigmap <- (1+xi[1])^((xi[1]+1)/2)
			rnormp(n,mu=0,sigmap=sigmap,p=p)
		}
	}	
	if(family=="Sinh-normal"){
	    if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
		rvg <- function(n){
		      asinh(rnorm(n)*xi[1]/2)
		}
	}	
	
	if(family=="Sinh-t"){
	    if(xi[1]<=0 | xi[2]<=0) stop("the extra parameters must be positive!!",call.=FALSE)
		rvg <- function(n){
		      asinh(rt(n,df=xi[2])*xi[1]/2)
		}
	}
	
	if(family=="Hyperbolic"){
	    if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
		rvg <- function(n){
               sqrt(rgig(n, lambda=1, chi=1, psi=xi[1]*xi[1]))*rnorm(n)
		}
	}
	
	if(family=="Slash"){
	    if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
		rvg <- function(n){
               sqrt(1/rbeta(n,shape1=xi[1],shape2=1))*rnorm(n)
		}
	}
	rvg(n)
}
