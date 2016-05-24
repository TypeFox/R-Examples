##
## Run DPM model
##

dpmixsim <-
function(x, M=1, a=1, b=1, upalpha=1, a0=2, b0=2, maxiter=4000, rec=3000,  fsave=NA, kmax=30, nclinit=NA, minvar=0.001)
{
    ## run dpmodel
    runif(1) # just for generating the seed if it doesn't exit
    n    <- length(x)
    cat("simulation length:",n,"\n")
    if(is.na(nclinit)) {
      sinit <- seq(0,n-1)
	    njinit <- rep(1,n)
  	  njlen <- n
    }
    else {
      ## initial custom initialization of clusters using cluster::clara
      sx <- clara(x, nclinit)
      sinit <- sx$clustering - 1 
	    njinit <- table(sinit)[]
  	  njlen <- length(njinit)
    }
    ## 
    ptm  <- proc.time()
    res  <- .C("gibbsdpm",
			  as.double(x),
	  		as.integer(sinit),
	  		as.integer(njinit),
	  		as.integer(njlen),
	  		as.integer(n),
	  		as.double(M),
        as.double(a),
        as.double(b),
        as.double(a0),
        as.double(b0),
        as.double(minvar),
        as.integer(upalpha),
        as.integer(maxiter),
        as.integer(rec),
        as.integer(kmax),
        krec   = as.integer(rep(0,rec)), # number of simulated components
        wrec   = as.double (rep(0,rec*kmax)), # weights
        phirec = as.double (rep(0,rec*kmax)),
        varrec = as.double (rep(0,rec*kmax)))
    cat("\ntime of gibbdpm: ", (proc.time() - ptm)[1]/60,"\n")
    res <- list(krec=res$krec,wrec=res$wrec,phirec=res$phirec,varrec=res$varrec)
    if(!is.na(fsave)) {
      cat("saving simulation ", fsave, "...")
      save(res, file = fsave)
      cat("\n")
    }
    invisible(res)
}

