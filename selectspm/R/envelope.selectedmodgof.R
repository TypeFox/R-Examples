envelope.selectedmodgof<- function(Y, fun=NULL, nrank=1,nsim=99,dimyx=c(128,128),...){

        x<- Y # the argument should be named Y in order to regeister as an envelope S3 method
	
	simu.model <- x$best.model
	cual <- class(simu.model)
	
	 bw <- x$best.sigma
	 
	 # envueltas para un HPC 
	 if("ecespa.minconfit"%in%cual & !is.na(bw) ){
	   if(is.null(fun)) fun<- Kinhom
	   
	   lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
	   simu.model$lambda <- lambda
	   result <- envelope(x$pp, fun, sigma=bw,
	        simulate=expression(rIPCP(simu.model)),
		savefuns =TRUE,nrank=nrank, nsim=nsim,...)
	}
	 # envueltas para un PC
        if("ecespa.minconfit"%in%cual & is.na(bw) ){
	   if(is.null(fun)) fun<- Kest
	  lambda <- predict(ppm(x$pp), type = "trend")
	  simu.model$lambda <- lambda
	  result <- envelope(x$pp, fun,
	        simulate=expression(rIPCP(simu.model)),
		savefuns =TRUE,nrank=nrank, nsim=nsim,...)
		
	}
	 # envueltas para un HPP
        if("im"%in%cual  & !is.na(bw) ){
	   if(is.null(fun)) fun<- Kinhom
	  
	   lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
  	   result <- envelope(x$pp, fun, sigma=bw,
  		  simulate=expression(rpoispp(lambda)),
		  savefuns =TRUE,nrank=nrank, nsim=nsim,...)
	}
	if("im"%in%cual & is.na(bw)){
	   if(is.null(fun)) fun<- Kest
	   result <- envelope(x$pp, fun,
	   savefuns =TRUE,nrank=nrank, nsim=nsim,...)
	}
	
	return(result)
}