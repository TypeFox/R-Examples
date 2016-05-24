"rstats" <-
function(RSobj,userfunc,...)
{
    obj.name <- deparse(substitute(RSobj))
    if (!(class(RSobj)=="RSmpl" || class(RSobj)=="RSmplext")){
         err.text<-paste(obj.name," is not a sample object - see help(\"rsextrobj\")",sep ="",collapse="")
         stop(err.text)
    }

    # extracts simulated matrices into three dimensional array sim
    n_tot  <- RSobj$n_tot
    n      <- RSobj$n
    k      <- RSobj$k
    nwords <- c(trunc((k+31)/32))

    # store coded simulated matrices in list with n_eff+1 elements
    sim<-split(RSobj$outvec,gl(n_tot,n*nwords))


    # decode simulated matrices and apply user function
    #RET<-unlist(lapply(sim,rsunpack,n,k,nwords,userfunc))
    RET<-lapply(sim,rsunpack,n,k,nwords,userfunc,...)
    RET
}

