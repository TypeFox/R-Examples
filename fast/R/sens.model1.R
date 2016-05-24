example_model1 <- function(par,a, output=c("model", "analytical sensitivities")){
	#Saltelli and Sobol 1995: About the use of rank transformation...
        #Davis and Rabinowith 1984
	#Example: example_model1(par=c(0.5,0.5,0.5),a=c(1,1,1))

	output <- match.arg(output, several.ok=FALSE)
	stopifnot(mode(a)=="numeric")
	stopifnot(all(a>=0))

	if(output=="model"){
		stopifnot(all(par<=1 & par>=0))
		stopifnot(length(par)==length(a))
		g <- function(a,par){
			#to.ret <- outer(X=par,Y=a, FUN=function(par,a){
			#           (abs(4*par-2)+a)/(1+a)
			#})
			return((abs(4*par-2)+a)/(1+a))
		}
		#return(apply(g(a,par), FUN=prod, MARGIN=1))
		return(prod(g(a,par)))
	} else {
                partial.var <- 1/(3*(1+a)^2)
		total.var <- prod(1+partial.var)-1
		return(partial.var/total.var)
	}


}
