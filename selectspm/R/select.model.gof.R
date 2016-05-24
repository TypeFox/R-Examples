
# Es el mismo codigo (expurgado) de "select.model.gofEXTRA_2.R"
# Permite suministrar parametros alternativos para la optimizacion
# compara los ajustes de varios metodos de optimizacion

select.model.gof <-
function (pp, sigmas, r, nlarge = 10000, q = 1/4, p = 2, correction = "trans",
              sigma2=NULL, rho=NULL, lower=NULL, upper=NULL, parscale=c(1,1),
	      dimyx=c(128,128), nsim=99, seed=1, correct.lambda=10) 
{
    len.sig <- length(sigmas)
    pprho <- pp$n/area.owin(pp$window)

    if (is.null(lower)) lower <- c(sigma2 = r[2]/10, rho = 1/area.owin(pp$window))
    if (is.null(upper)) upper <- c(sigma2 = (max(r) * 4)^2, rho = pprho)
    if(!is.null(sigma2)) sigma2.0 <- sigma2 else sigma2.0 <- (upper["sigma2"] - lower["sigma2"])/2
    if(!is.null(rho)) rho.0 <- rho  else rho.0 <- (upper["rho"] - lower["rho"])/2
    
    # alternative parscale
    parscale2<-c(max(upper),min(lower))
    
    
    
   # lists to store inhomogeneous (*H*eterogeneous) Poisson and Poisson cluster models
   HPPs = list()
   HPCs = list()
   models <- list() # List to store models for subsequent simulation
    
    HPC.gofs<-list()
    PC.gofs<-list()
    
    for (i in 1:len.sig) {
        cat(paste("evaluating bandwith", i, "of ",len.sig,"\n\n"))
    
        lambda <- density.ppp(pp, sigma = sigmas[i], at = "points")
	
	 
	 # correction of lambda values ==0.
         # We assign a fraction (1/correct.lambda) of the smalest lambda value.
	 # This error will occassionaly occur for small bandwiths in ppp with few points
	 lambda[lambda==0] <- min(lambda[lambda!=0])/correct.lambda
	 
	#realiza tres tipos de optimizacion: dos BFGS y una NElder-MEad
        hpc.model <- list(
	    hpc.model1 <- ipc.estK2(pp, lambda = lambda, correction = correction, 
             r = r, nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, 
             rho = rho.0, method = "L-BFGS-B", lower = lower, 
            upper = upper, control = list(parscale = parscale)),
	    
           hpc.model2 <- ipc.estK2(pp, lambda = lambda, correction = correction, 
             r = r, nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, 
            rho = rho.0, method = "L-BFGS-B", lower = lower, 
            upper = upper, control = list(parscale = parscale2)),
	    
	  hpc.model4 <- ipc.estK2(pp, lambda = lambda, correction = correction, 
             r = r, nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, 
             rho = rho.0)
	)
	
	
	 # recalculate  lambda at pixels in order to simulate HPPs and  HPCs 
	 lambda <- density.ppp(pp, sigma = sigmas[i], at = "pixels", dimyx=dimyx)
	
	 # check that the fitted parameters are within the proposed range 
	 hpc.par.ok<- sapply(hpc.model, function(x){
	        
	        x$sigma2>=lower[1] & x$sigma2 <=upper[1] & x$rho>=lower[2] & x$rho <=upper[2]
	        })
	 hpc.model.ok <- hpc.model[hpc.par.ok]
	 

	 
	 
	 # Among the fitted models with parameters within the allowed range, 
	 #chose the one whose simulations produce the best fit (the smalest "u" value of the GoF test) 
	
	 if (length(hpc.model.ok)>0){ 
	     test.HPC<- list()
	     gof.test.HPC<-list()
	     for (m in 1:length(hpc.model.ok)){
	         cat(paste("evaluating adjustment for HPC model", m, "of ",length(hpc.model.ok),"\n"))
	     
	         hpc.model.ok[[m]]$lambda <- lambda # plug lambda into the model in order to simulate them with rIPCP
		 set.seed(seed)
	 	 test.HPC[[m]] <-envelope(pp, Kinhom, sigma = sigmas[i], r=r, nsim=nsim, 
	                                     simulate=expression(rIPCP(hpc.model.ok[[m]])),savefuns=TRUE,
					     correction=correction, nlarge=nlarge)
		gof.test.HPC[[m]]<-LF.gof(test.HPC[[m]])
		
		}
	
	        # evaluate the goodnes of fit of every HPC
	       gof.test.HPC.u<- sapply(gof.test.HPC, function(x)x$u)
	       
	       # Store the best model and all the related stuff
	       HPCs[[i]] <- test.HPC[[which.min(gof.test.HPC.u)]]
	       hpc.model.ok <-hpc.model.ok[[which.min(gof.test.HPC.u)]]
	         
	       HPC.gofs[[i]] <-gof.test.HPC[[which.min(gof.test.HPC.u)]]
	
	 }
	 
	  if (length(hpc.model.ok)<1) {
         	HPCs[[i]] <-NA
		HPC.gofs[[i]] <- list(u=NA, p=NA, na.count.by.r=NA)
		hpc.model.ok <- NA
	  }
	
         
	 cat(paste("evaluating adjustment for HPP model","\n"))
	# store HPP and the best  HPC 
	set.seed(seed)
	HPPs[[i]] <- envelope(pp, Kinhom, sigma = sigmas[i], r=r, nsim=nsim, 
	                                   simulate=expression(rpoispp(lambda)),savefuns=TRUE,
					   correction=correction, nlarge=nlarge)
        models[[i]] <- lambda
        models[[i + len.sig]] <- hpc.model.ok
       
    }
    
    # Chose the best homogeneous PC  model among three  implementations of the optimization algorithm
    pc.model <-list(
        pc.model1<- ipc.estK2(pp, correction = correction, r = r, 
            nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, rho = rho.0, 
            method = "L-BFGS-B", lower = lower, upper = upper, control = list(parscale = parscale)),
	pc.model2<- ipc.estK2(pp, correction = correction, r = r, 
            nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, rho = rho.0, 
            method = "L-BFGS-B", lower = lower, upper = upper, control = list(parscale = parscale2)),
	
	pc.model4<- ipc.estK2(pp, correction = correction, r = r, 
            nlarge = nlarge, p = p, q = q, sigma2 = sigma2.0, rho = rho.0)
	)
	
	
	 # check that the fitted parameters are within the proposed range 
	 pc.par.ok<- sapply(pc.model, function(x){
	        
	        x$sigma2>=lower[1] & x$sigma2 <=upper[1] & x$rho>=lower[2] & x$rho <=upper[2]
	        })
	 pc.model.ok <- pc.model[pc.par.ok]
	 
	 
	 # Among the fitted models with parameters within the allowed range, 
	 #chose the one whose simulations produce the best fit (the smalest "u" value of the GoF test) 
	 
	#  homogeneous lambda required to simulate with  rIPCP and rpoispp
	 lambda.homo <- predict(ppm(pp), type = "trend", locations=lambda)  
	
	if (length(pc.model.ok)<1) {
         	PCP <-NA
		PC.gof<- list(u=NA, p=NA, na.count.by.r=NA)
		pc.model.ok <- NA
	  }
	 if (length(pc.model.ok)>0){ 
	     test.PC<- list()
	     gof.test.PC<-list()
	     
	     
	     for (m in 1:length(pc.model.ok)){
	     cat(paste("evaluating adjustment for PC model", m, "of ",length(pc.model.ok),"\n"))
	        pc.model.ok[[m]]$lambda <-lambda.homo # la metemos en el modelo para poder simular con rIPCP
		set.seed(seed)
	 	 test.PC[[m]] <-envelope(pp, Kinhom, sigma = sigmas[i], r=r, nsim=nsim, 
	                                     simulate=expression(rIPCP(pc.model.ok[[m]])),savefuns=TRUE,
					     correction=correction, nlarge=nlarge)
		gof.test.PC[[m]]<-LF.gof(test.PC[[m]])
		
		}
	
	        # evaluathe the GoF of every PC
	       gof.test.PC.u<- sapply(gof.test.PC, function(x)x$u)
	       
	      # store the best model and all the related stuff
	       PCP <- test.PC[[which.min(gof.test.PC.u)]]
	       pc.model.ok <-pc.model.ok[[which.min(gof.test.PC.u)]]
	          
	       PCP.gof<-gof.test.PC[[which.min(gof.test.PC.u)]]
	
	 }
	
          	  
	 	
         models[[(2 * len.sig) + 1]] <- pc.model.ok
         
    
        # FINALY, "adjust " the homogeneopus Poisson model etc
	cat(paste("evaluating adjustment for P model","\n"))
	set.seed(seed)
	 P <- envelope(pp, Kest,  r=r, nsim=nsim, simulate=expression(rpoispp(lambda.homo)),
	             savefuns=TRUE, correction=correction, nlarge=nlarge)
    
   
    models[[(2 * len.sig) + 2]] <- lambda.homo
    
    
    # compute and/or store gof's
    gofs <- list()
    # first, for HPP's
    for ( i in 1:length(sigmas)){
           gofs[[i]] <-LF.gof( HPPs[[i]])
    }
   #then for  HPC's
   for ( i in 1:length(sigmas)){
	      gofs[[length(sigmas)+i]] <-HPC.gofs[[i]] # previously calculated in the evaluation of the best fiting optimization
	      
	   }
   #then for the PC
   gofs[[(2*length(sigmas))+1]] <- PCP.gof
      
    #finaly for the  P
   gofs[[(2*length(sigmas))+2]] <- LF.gof(P)
   

   # selecct and return the best model and all related stuff
   gof.u<- sapply(gofs, function(x)x$u)
   nombres.modelos<- c(paste("HPP_sg_", sigmas), paste("HPC_sg_", sigmas), "PC", "P")
   
   names(gof.u) <- nombres.modelos
   best.gof<- which.min(gof.u)
   best.sigma<-NA # to store best bandwith for simulating from the fitted model
   if(best.gof<=length(sigmas)){
               best.envelopes<-HPPs[[best.gof]]
	       best.sigma<-sigmas[best.gof]
       }
   if(best.gof>length(sigmas) & best.gof<=(2*length(sigmas))){
                best.envelopes<-HPCs[[best.gof-len.sig]]
                best.sigma <- sigmas[best.gof-len.sig]
	}
   if(best.gof>(2*length(sigmas)) & best.gof<((2*length(sigmas)) +2)   ) best.envelopes<-PCP
   if(best.gof>((2*length(sigmas)) +1)) best.envelopes<-P
   result<- list(gof.u=gof.u, best.gof =gof.u[best.gof], best.model=models[[best.gof]], gof=gofs[[best.gof]], envelopes=best.envelopes, pp=pp,best.sigma=best.sigma)
   class(result) <- c("selectedmodgof",class(result))
   
   return(result)
   }
   
    