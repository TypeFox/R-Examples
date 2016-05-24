print.thetaMSAR.VM <-
function(x,digits=4,...) {
    if (!is.null(attr(x,'n_par'))) {
        n_par <- attr(x,'n_par')
    } else {
        n_par <- 0
    }

    thetadim <- dim.thetaMSAR(x)
    order <- attributes(x)$order
    label=attributes(x)$label
    M=attributes(x)$NbRegimes
    d=attributes(x)$NbComp
    if(label=='NULL'){label='HH'}
    ncov.trans = length(x$par.trans[1,])-1
    ncov.emis = length(c(x$par.emis[[1]]))
    if (thetadim[2]>=2) {
        x <- as.thetaMSAR(x,label=label,ncov.emis=ncov.emis,ncov.trans=ncov.trans)  # this adds dimnames if they don't exist
        sel <- 1:n_par
    } else {
        sel <- TRUE
    }
    class(x) <- NULL
    dimname <- dimnames(x)
    if(order>0){cat("AUTOREGRESSIVE COEFFICIENTS (kappa)\n")
    	if(d==1){
    		print(x$kappa,digits)
    		cat("\n")
    	} else {
    		for(i in 1:M){
    			cat(paste("Regime",i,"\n"))
    			for(j in 1:(order+1)){
    				cat(paste("kappa",j,sep=""))
	    			print(x$kappa[i,j],digits)
	    			cat("\n")
    			}
    		}
    	}
    }
    if (order==0) {
    	cat("MEAN\n") 
    	print(x$mu,digits)
    } else {
    	cat("Intercept\n")
    	print(x$mu,digits)
   
#    cat("MEAN\n")
#		A = NULL
#    	for(m in 1:M){ 
#    		for (o in 1:order) {
#    			A[m,] = A[[m]]+x$A[m,o]
#    		}
#    	}
#        print(solve(diag(d)-A[[m]])%*%x$mu,digits)
	}
    cat("\nINITIAL DISTRIBUTION OF THE MARKOV CHAIN")
    print(x$prior,digits)
    cat("\nTRANSITION PROBABILITY MATRIX\n")
    print(x$transmat,digits)
    if(label=="HN" || label=="NN"){
    	cat('\nEMISSION PARAMETERS\n')
    	for(i in 1:M){
    		cat("Regime",i,'\n')
    		print(x$par.emis[[i]],digits)
    		cat('\n')
    	}
    }
    if(label=='NH' || label=='NN'){
    	cat("\nTRANSITION PARAMETERS\n")
    	print(x$par.trans,digits)
    }
}
