"getClpConstr" <-
function (clp, clp0, clpequspec, ncomp, compnames) {

    ## doing it this way because for mass spectra case 
    ## clp0 >> ncomp 
    ## and often length(clpequspec) < ncomp
    clp0mat <- matrix(0, length(clp), ncomp, dimnames = list(c(),compnames))
    clpRem <- matrix(0, length(clp), length(clpequspec))
    clpMod <- matrix(0, length(clp), length(clpequspec))

    if (length(clp0) != 0) 
        for (i in 1:length(clp0)){
            to0 <- intersect(which(clp >= clp0[[i]]$low), which(clp <= 
                  clp0[[i]]$high))
	    if(length(to0>0)) {
	       for(j in 1:length(to0)) 
            	clp0mat[to0[j], clp0[[i]]$comp] <- 1
            }
	}
    if (length(clpequspec) != 0) 
        for (i in 1:length(clpequspec)) {
            clptoequ <-  intersect(which(clp >= clpequspec[[i]]$low), 
                  which(clp <= clpequspec[[i]]$high))
           
            clpRem[clptoequ, i] <- clpequspec[[i]]$to
            clpMod[clptoequ, i] <- clpequspec[[i]]$from
         }
    list(clp0mat = clp0mat, clpRem = clpRem, clpMod = clpMod)

}
