flaremixEM<-function(y, x, lambda = NULL, beta = NULL, sigma = NULL, alpha = NULL, nu=NULL,
		epsilon = 1e-04, maxit = 10000, verb = FALSE, restart=50){
j=1
while(j<=length(nu)){
temp=try(try.flare(y=y, x=x, lambda = lambda, beta = beta, sigma = sigma, alpha = alpha, nu=nu[j],
		epsilon = epsilon, maxit = maxit, verb = verb, restart=restart),silent=TRUE)

if(class(temp)=="try-error") j=j+1 else j=2^100
}



if(j==(length(nu)+1)) stop(paste("Too many attempts. Select a different barrier constant.","\n"))


temp

}

