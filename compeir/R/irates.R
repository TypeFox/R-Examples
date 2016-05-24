irates <-
function(	data, 
					time.code=NULL, 
					no.event.code=NULL,
					no.event.lab = NULL, 
					event.lab=NULL, 
					covar.lab=NULL,
					full.sample.lab = "Full sample",
					ci.level = .95,
					ci.fun = "log"){

		quant <- 1-(1-ci.level)/2

		if(!(ci.fun %in% c("lin", "log"))) stop(paste("ci.fun", ci.fun, "is not implemented in irates"))
		
		if(is.null(time.code)) time.code = colnames(data)[1]
		if(is.null(no.event.code)) no.event.code = colnames(data)[2] 
		
		time.code = as.character(time.code)
		no.event.code = as.character(no.event.code)
		
		covar.code = rownames(data)
		full.sample.code = "full.sample"
		
		fullsample = c(rownames(data), full.sample.code)
		event.code = colnames(data)[-c(which(colnames(data) == time.code), which(colnames(data) == no.event.code))]
		
		if(!is.null(event.lab) && length(event.lab) != length(event.code)){
			stop("Length event.lab != length event")
			}
		if(!is.null(covar.lab) && length(covar.lab) != length(rownames(data))){
			stop("Length covar.lab != length cov.code")
			}

		empty.data <- as.data.frame(matrix(NA, nrow=length(fullsample), ncol=length(event.code)))
		names(empty.data) <- event.code
		rownames(empty.data) <- fullsample
		ir <- var <- conf.upper <- conf.lower <- empty.data
    
	    n <- data[fullsample, event.code]
		ir <- n/data[fullsample, time.code]

		N = rowSums(data[-which(names(data) == time.code)])
		
		if(any(fullsample == full.sample.code)){
			if(length(event.code) == 1){
				if(length(fullsample) == 1){
					n <- sum(data[,event.code])
					ir <- sum(data[,event.code])/sum(data[,time.code])
					} else{
						n[which(is.na(n))] <- sum(data[,event.code])
						ir[which(is.na(ir))] <- sum(data[,event.code])/sum(data[,time.code])
						}
			} else{	
					n[which(rownames(n) == "NA"), ] <- colSums(data[,event.code])
					ir[which(rownames(ir) == "NA"), ] <- colSums(data[,event.code])/sum(data[,time.code])
					}
			}

	   	var <- ir^2/n

		   
   	    if ( ci.fun == "log" )
      	{ 
        	conf.lower <- ir*exp(-qnorm(quant)/sqrt(n))
        	conf.upper <- ir*exp(+qnorm(quant)/sqrt(n))
	    }
      	if ( ci.fun == "lin" )
      	{
        	conf.lower <- ir-qnorm(quant)*sqrt(var)
        	conf.upper <- ir+qnorm(quant)*sqrt(var)
        }
 
 		### labels
 		no.event.name = no.event.code
 		if(!is.null(no.event.lab)) no.event.name = no.event.lab  		
 		event.name = event.code
 		if(!is.null(event.lab)) event.name = event.lab
  		
  		covar.name = covar.code
  		if(!is.null(covar.lab)) covar.name = covar.lab
  		
  		fullsample.name = full.sample.code
  		if(!is.null(full.sample.lab)) fullsample.name = full.sample.lab
  		  		
		### final data.frames
		names(N) = covar.code
  		n = as.data.frame(n, row.names=fullsample)
  		names(n) = event.code
  		ir = as.data.frame(ir, row.names=fullsample)
  		names(ir) = event.code		
  		var = as.data.frame(var, row.names=fullsample)
  		names(var) = event.code
  		conf.lower = conf.lower=as.data.frame(conf.lower, row.names=fullsample)
  		names(conf.lower) = event.code
  		conf.upper = as.data.frame(conf.upper, row.names=fullsample)
  		names(conf.upper) = event.code
  
  		res <- list(
    	          	ir = ir,
	              	var=var,
        	      	conf.lower=conf.lower,
            	  	conf.upper=conf.upper,
              		no.event.code = no.event.code,
              		no.event.lab = no.event.name,
              		event.code = event.code,
              		event.lab = event.name,
              		covar.code = covar.code,
              		covar.lab = covar.name,
              		full.sample.code = full.sample.code,
              		full.sample.lab = fullsample.name,
              		N=N,
              		n=n,
              		ci.level=ci.level,
              		ci.fun=ci.fun)
              
  class(res) <- "irates"

  return(res)
}

