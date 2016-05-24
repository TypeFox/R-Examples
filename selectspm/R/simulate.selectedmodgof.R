simulate.selectedmodgof<- function(object, nsim=99, seed=1,dimyx=c(128,128),...){

	x<-object # hay que usar "object" as argument in order to register as a S3 simulate method

        result<-list()
	simu.model <- x$best.model
	cual <- class(simu.model)
	
	 bw <- x$best.sigma
	 
	 # envueltas para un HPC 
	 if("ecespa.minconfit"%in%cual & !is.na(bw) ){	   
	   lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
	   simu.model$lambda <- lambda
	    for( i  in 1:nsim){
	           set.seed(seed)
	           progressreport(i,nsim)
	   	   result[[i]] <- rIPCP(simu.model)
	   }
	}
	
	 # envueltas para un PC
        if("ecespa.minconfit"%in%cual & is.na(bw) ){
	  lambda <- predict(ppm(x$pp), type = "trend")
	  simu.model$lambda <- lambda
	 for( i  in 1:nsim){
	           set.seed(seed)
	           progressreport(i,nsim)
	   	   result[[i]] <- rIPCP(simu.model)
	   }
	}
	
	 # envueltas para un HPP
        if("im"%in%cual  & !is.na(bw) ){
	 lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
  	 for( i  in 1:nsim){
	           set.seed(seed)
	           progressreport(i,nsim)
  	           result[[i]] <- rpoispp(lambda)
	   }
	}
	
	# envueltas para un P
	if("im"%in%cual & is.na(bw)){
	   lambda <- intensity(x$pp)
	   ventana<- x$pp$window
	   for( i  in 1:nsim){
	           set.seed(seed)
	           progressreport(i,nsim)
	           result[[i]] <- rpoispp(lambda, win=ventana)
	   }
	}
	
	return(result)
}