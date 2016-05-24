getPTLparams <- function(x1,x2){

	#cat("fitting PTL parameters...")
        if(!length(x1)==length(x2)) stop("vectors supplied to getPTLparams() must be of equal length")
        bin_limits <-  seq(from=-2,to=2,length=51)

        d <- sqrt(sum((x1-x2)^2,na.rm=TRUE))
        vars <- x1-x2
        vars <- vars[!is.na(vars)]
	#cat(vars)
	#cat("\n")
	#cat(paste("d:",d))
        beta_nls <- NA
        alpha_nls <- NA
        gamma_nls <- NA

        if((length(vars)/length(x1))>0.1){
                beta_hat <- estimateB(vars)

                observed_table <- table(cut(vars,breaks=bin_limits))
                observed_counts <- rep(NA,length(bin_limits)-1)
                for(i in 1:length(observed_counts)){
                        observed_counts[i] <- observed_table[[i]]
                }

		nls_fit <- NULL
	
		#nls.error <- function(e){
		#	cat(paste(e,"triggered by d=",d,"and beta:",beta_hat,"\n"))
		#	cat(paste("length(observed_counts)=",length(observed_counts),"\n"))
		#	#e
		#	NULL
		#}

		PTL.nls <- function(alpha,beta,gamma){getPTLExpectedCounts(alpha=alpha,beta=beta,gamma=gamma,bin_limits=bin_limits,ntrials=length(vars))}

                #nls_fit <- tryCatch(nls_fit <- nls(observed_counts ~ getPTLExpectedCounts(alpha=linear_param,beta=scale_param,gamma=poly_param,bin_limits=bin_limits,ntrials=length(vars)),start=list(scale_param=beta_hat/2,linear_param=0,poly_param=0),control=nls.control(maxiter=100,warnOnly=TRUE)),error=nls.error)
		nls_fit <- tryCatch(nls(observed_counts ~ PTL.nls(alpha=linear_param,beta=scale_param,gamma=poly_param),start=list(scale_param=beta_hat/2,linear_param=0,poly_param=0)),error=function(e) NULL)
		#nls_fit <- nls(observed_counts ~ getPTLExpectedCounts(alpha=linear_param,beta=scale_param,gamma=poly_param,bin_limits=bin_limits,ntrials=length(vars)),start=list(scale_param=beta_hat/2,linear_param=0,poly_param=0),control=nls.control(maxiter=100,warnOnly=TRUE))
			
        
                if(!is.null(nls_fit)){
			#if(!nls_fit$convInfo$isConv) cat("nls fit failed \n")
	                if(nls_fit$convInfo$isConv){
	                        #cat("PTAL model fitted successfully \n ")
	                        beta_nls <- coef(nls_fit)[1]
	                        alpha_nls <- coef(nls_fit)[2]
	                        gamma_nls <- coef(nls_fit)[3]
	                }
 		}
	}

        list(dist=d,beta=beta_nls,alpha=alpha_nls,gamma=gamma_nls)

}

