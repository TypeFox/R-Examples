correction.redundancy <-function(r,HXmax,Herror,finite){
	
    
    out<-.C("correctionredundancy",
            HXmax=as.double(HXmax),
            Herror=as.double(Herror[finite]),
            r=as.double(r),
            lengthRedundancy=as.double(length(r)))
    
	return(out$r)
	}

