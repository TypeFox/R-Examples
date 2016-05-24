# TODO: Add comment
# 
# Author: Giorgio Spedicato
###############################################################################


summary.cwrObj <-function(object,...)  #metodo summary per cwr objects
{
	cat("cwrObj : ", deparse(substitute(object)), "\n")        #stampa il nimer
	cat("-------------","\n", "\n")
	cat("X component","\n")    #per ogni gruppo i momenti stimati
	cat("muX","\n")
	print(object$muX)
	cat("sigmaX", "\n")
	print(object$SigmaX)
	cat("\n","Y component","\n")
	cat("muY","\n")
	print(object$muY)
	cat("sigmaY", "\n")
	print(object$SigmaY)
	cat("\n", "Weights (beta)", "\n")
	print(object$weightsY)
	#fa i parametri delle rette rapperesntare parametri rette
	if((size(object$muX,1)==1)&(size(object$muY,1)==1))
	{ 
		cat("---------------------", "\n")
		cat("Regression equations for each cluster","\n")
		for(i in 1:(size(object$muX, 2)))
		{
			cat("Equation for group ", i, " is: y=", object$muY[,i],"+",object$weightsY[1,1,i],"x","\n")
		}
	}
}

#metodo print x oggetti CWM
print.cwrObj<-function(x,...)   #metodo print per oggetti cwr
{
	cat("cwrObj : ", deparse(substitute(x)), "\n")     #stampa il nome
	cat("-------------","\n", "\n")
	cat("Number of groups: ", size(x$muX, 2), "\n")      #numero e pesi (priors ) dei gruppi
	cat("Relative weights (priors): ", x$priorC, "\n")
	cat("-------------","\n", "\n")
	cat("logLik : ", x$logLik, "\n")  #logLik aik e bic
	cat("AIC : ", x$aic, "\n")
	cat("BIC : ", x$bic, "\n")
	cat("-------------","\n", "\n")
	cat("computation time : ", x$timeTotal, "\n")
}

#plot method

plot.cwrObj<-function(x,...)    #metodo plot per cwr objects . funziona solo se stanno in R^2
{
	if(!((size(x$muX,1)==1)&(size(x$muY,1)==1))) stop("Error! Input and output shall be unidimensional")
	with(x, plot(Y~X, col=group, ...))
}

#logLik method

logLik.cwrObj<-function(object,...)
{	out<-numeric(1)
	out<-object$logLik
	return(out)
}

#predict method

predict.cwrObj<-function(object,...)
{
	group=apply(object$posteriors, 2, .maxIndex) #assigns group
	group=t(as.matrix(group)) #save group membership 
	return(as.numeric(group))
}


