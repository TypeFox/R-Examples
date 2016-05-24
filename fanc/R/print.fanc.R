print.fanc=function(x,digits = max(3, getOption("digits") - 3),num.result = 20,...){
		#check digits
		if(mode(digits)!="numeric")	stop('"digits" must be numeric.')
		if(length(digits) > 1)	stop('"digits" must be a scalar (1-dimensional vector).')
		if(as.integer(digits)!=digits)	stop('"digits" must be integer.')
		if(digits <= 0)	stop('"digits" must be positive integer.')

		#check num.result
		if(mode(num.result)!="numeric")	stop('"num.result" must be numeric.')
		if(length(num.result) > 1)	stop('"num.result" must be a scalar (1-dimensional vector).')
		if(as.integer(num.result)!=num.result)	stop('"num.result" must be integer.')
		if(num.result <= 0)	stop('"num.result" must be positive integer.')

		
		 rho <- x$rho
		 gamma <- x$gamma
		 colnames(rho) <- sapply(gamma,function(x) sprintf('%3.3f',x))
		 df <- x$df[,1]
		 rho0 <- cbind(df,rho)
		 colnames(rho0) <- c("df",colnames(rho))


#RESULT

    # cat("\nCall:", paste(deparse(x$call)), "\n\n")
	# cat("\nThe result of the AIC:\n")
    # cat("\nUniquenesses:\n"); print(diagPsiAIC,digits=digits);
    # print(LambdaAIC);
    # cat("\nrho:\n"); print(rhoAIC,digits=digits);
    # cat("\ngamma:\n"); print(gammaAIC,digits=digits);
    # cat("\n")
 	# invisible(x)
 	
 	cat("\nCall:", paste(deparse(x$call)), "\n")
    cat("\ngamma:\n"); print(gamma,digits=digits);
    cat("\nrho:\n"); print(rho,digits=digits);
    cat("\ncor.factor:", paste(x$cor.factor));
    cat("\n")
 	invisible(x)      
    
}
