envelope.selectedmod<- function(Y, fun=NULL, nrank=1,nsim=99,dimyx=c(128,128),...){
	
	x<- Y # the argument should be named Y in order to regeister as an envelope S3 method
	
	# que tipo de modelo es el mejor
	modelos<-c("P","HPC","PC","HPP")
	cual <-  modelos[which(!is.na(pmatch(modelos,names(x$best.dtheta))))]
	print(cual)
	sumario<-print(x)
         
	 simu.model <- x$best.model
	 # envueltas para un HPC
	 if(cual=="HPC"){
	   if(is.null(fun)) fun<- Kinhom
	   bw <-  sumario["HPC.bw"]
	   lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
	   simu.model$lambda <- lambda
	   result <- envelope(x$pp, fun, sigma=bw,
	        simulate=expression(rIPCP(simu.model)),
		savefuns =TRUE,nrank=nrank, nsim=nsim,...)
	}
	 # envueltas para un PC
        if(cual=="PC"){
	   if(is.null(fun)) fun<- Kest
	  lambda <- predict(ppm(x$pp), type = "trend")
	  simu.model$lambda <- lambda
	  result <- envelope(x$pp, fun,
	        simulate=expression(rIPCP(simu.model)),
		savefuns =TRUE,nrank=nrank, nsim=nsim,...)
		
	}
	 # envueltas para un HPP
        if(cual=="HPP"){
	   if(is.null(fun)) fun<- Kinhom
	   bw <- sumario["HPP.bw"]
	   lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
  	   result <- envelope(x$pp, fun, sigma=bw,
  		  simulate=expression(rpoispp(lambda)),
		  savefuns =TRUE,nrank=nrank, nsim=nsim,...)
	}
	if(cual=="P"){
	   if(is.null(fun)) fun<- Kest
	   result <- envelope(x$pp, fun,
	   savefuns =TRUE,nrank=nrank, nsim=nsim,...)
	}
	
	return(result)
}