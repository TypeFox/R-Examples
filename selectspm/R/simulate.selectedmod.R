simulate.selectedmod<- function(object,nsim=99, seed=1, dimyx=c(128,128),...){

        x<-object # hay que usar "object" as argument in order to register as a S3 simulate method

	# que tipo de modelo es el mejor
	modelos<-c("P","HPC","PC","HPP")
	cual <-  modelos[which(!is.na(pmatch(modelos,names(x$best.dtheta))))]
	print(cual)
	sumario<-print(x)
	
	result<-list()
         
	 simu.model <- x$best.model
	 # simulating  HPC
	 if(cual=="HPC"){
	   bw <-  sumario["HPC.bw"]
	   lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
	   simu.model$lambda <- lambda
	   for( i  in 1:nsim){
	           set.seed(seed)
	           progressreport(i,nsim)
	   	   result[[i]] <- rIPCP(simu.model)
	   }
	}
	 # simulating PC
        if(cual=="PC"){
	
	 lambda <- predict(ppm(x$pp), type = "trend")
	 simu.model$lambda <- lambda
	 for( i  in 1:nsim){
	           set.seed(seed)
	           progressreport(i,nsim)
	           result[[i]] <- rIPCP(simu.model)
	   }
	}
	 # simulating HPP
        if(cual=="HPP"){
	   bw <- sumario["HPP.bw"]
	   lambda <- density.ppp(x$pp, sigma=bw, dimyx=dimyx)
	   for( i  in 1:nsim){
	           set.seed(seed)
	           progressreport(i,nsim)
  	           result[[i]] <- rpoispp(lambda)
	   }
	}
	if(cual=="P"){
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
