# see vignette for notation

KFfit <- function(  param,              # model parameters
                    data,               # a matrix containing the data ie Y           
        		    ,...,               # delete this line and include other arguments to be passed
                                        #     to the function here
        		    optim=FALSE,        # set this to FALSE if not estimating parameters, otherwise,
        		                        #     if parameters have already been estimated, then set to TRUE
        		    history.means=FALSE,# whether to save a matrix of filtered E(Theta_t)
        		    history.vars=FALSE, # whether to save a list of filtered V(Theta_t)
        		    prior.mean=NULL,    # optional prior mean. column vector.    # Normally set 
        		    prior.var=NULL,     # optional prior mean. matrix.           # these inside the code
        		    fit=FALSE,          # whether to return a matrix of fitted values
        		    se.fit=FALSE,       # whether to return the standard error of the fitted values
        		    se.predict=FALSE,   # whether to return the prediction standard error =
        		                        #     se(fitted values) + observation variance V_t
        		    noisy=TRUE          # whether to print a progress bar, useful.
        		    na.rm=FALSE){       # whether to use NA handling. set to TRUE if any Y is NA 
		  
    start <- Sys.time()  
    
    if(se.predict & !se.fit){
        stop("Must have se.fit=TRUE in order to compute se.predict") # leads to a computational saving 
    }    
    
    #
    # Here, I would suggest creating dummy names for your paramters eg
    # sigma.obs <- exp(model.param[1])
    #
    #

    T <- dim(data)[1]   # ASSUMES OBSERVATIONS ARE IN A (T x m) matrix, ie row t contains the data for time t. 
                        # Note this is important for when KFadvance is called later     
    
    if(is.null(prior.mean)){
        Xpost <- DEFINE PRIOR MEAN HERE
    }
    else{
        Xpost <- prior.mean
    }
    
    if(is.null(prior.var)){   
        Vpost <- DEFINE PRIOR VARIANCE HERE
    }
    else{
        Vpost <- prior.var
    }    
    
    if (history.means){
        Xrec <- Xpost
    }
    if(history.vars){
        Vrec <- list()
        Vrec[[1]] <- Vpost
    }
    
    # delete or complete the following rows as necessary, also appears in the loop that follows
    A <- IF A IS FIXED OVER TIME THEN DEFINE IT HERE
    B <- IF B IS FIXED OVER TIME THEN DEFINE IT HERE
    C <- IF C IS FIXED OVER TIME THEN DEFINE IT HERE
    W <- IF W IS FIXED OVER TIME THEN DEFINE IT HERE   
    D <- IF D IS FIXED OVER TIME THEN DEFINE IT HERE
    E <- IF E IS FIXED OVER TIME THEN DEFINE IT HERE
    F <- IF F IS FIXED OVER TIME THEN DEFINE IT HERE
    V <- IF V IS FIXED OVER TIME THEN DEFINE IT HERE
    
    loglik <- 0 
    fitmat <- c() 
    sefitmat <- c()
    sepredictmat <- c()
    
    if(noisy){
        pb <- txtProgressBar(min=1,max=T,style=3)
    }
      
    for(t in 1:T){ 
    
        # delete or complete the following rows as necessary
        A <- IF A IS TIME-VARYING THEN DEFINE IT HERE
        B <- IF B IS TIME-VARYING THEN DEFINE IT HERE
        C <- IF C IS TIME-VARYING THEN DEFINE IT HERE
        W <- IF W IS TIME-VARYING THEN DEFINE IT HERE   
        D <- IF D IS TIME-VARYING THEN DEFINE IT HERE
        E <- IF E IS TIME-VARYING THEN DEFINE IT HERE
        F <- IF F IS TIME-VARYING THEN DEFINE IT HERE
        V <- IF V IS TIME-VARYING THEN DEFINE IT HERE
        
        # this bit calls KF advance      
        new <- KFadvance(obs=data[t,],oldmean=Xpost,oldvar=Vpost,A=A,B=B,C=C,D=Dt,E=E,F=F,W=W,V=V,marglik=TRUE,log=TRUE,na.rm=na.rm)
        
        Xpost <- new$mean
		Vpost <- new$var
		
		if(t==1){ # used when this function is called iteratively one step at a time
		    running.mean <- Xpost
		    running.var <- Vpost
		}
		
		loglik <- loglik + new$mlik
			
		
		if (history.means){
            Xrec <- cbind(Xrec,Xpost)
        }
        if(history.vars){
            Vrec[[t+1]] <- Vpost # since first entry is the prior
        }
        
        if(fit){
            fitmat <- cbind(fitmat,D%*%Xpost + E)
        }
        if(se.fit){
            sefitmat <- cbind(sefitmat,sqrt(diag(D%*%Vpost%*%t(D))))
        }
        if(se.predict){
            sepredictmat <- cbind(sepredictmat,sqrt((sefitmat[,ncol(sefitmat)])^2+diag(F%*%V%*%t(F))))
        }
        
        if(noisy){
            setTxtProgressBar(pb,t) 
        }
    }
    
    if(noisy){        
        close(pb)
    }
    end <- Sys.time()
    
    if(noisy){        
        cat("Time taken:",difftime(end,start,units="secs")," seconds.\n")
    }
    
    if(optim){
        return(-loglik) # just return the -log likelihood if in parameter estimation mode
    }
    else{
        retlist <- list(mean=Xpost,var=Vpost,mlik=loglik,data=data,running.mean=running.mean,running.var=running.var)
        if(history.means){
            retlist$history.means <- Xrec
        }
        if(history.vars){
            retlist$history.vars <- Vrec
        } 
        if(fit){
            retlist$fit <- fitmat
        }
        if(se.fit){
            retlist$se.fit <- sefitmat
        }
        if(se.predict){
            retlist$se.predict <- sepredictmat        
        }
        return(retlist)
    }
}

