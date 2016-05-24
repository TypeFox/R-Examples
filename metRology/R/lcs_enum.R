#LCS by full enumeration
#... probably stupid
#returns an array or matrix, each row of which contains the indices of a 
#consistent subset of smallest size
LCS<-function(x,u,p=0.05, method="enum", simplify=FALSE, verbose=FALSE) {
	
	if (method != "enum") 
        	warning(gettextf("method = '%s' is not supported. Using 'enum'", 
        	method), domain = NA)

	#Weighted mean internal function
	wtm<-function(x,u) sum(x/u^2)/sum(1/u^2)

	#Consistency check function:
	consistent<-function(x,u,theta,p) {
		pchisq(sum(((x-theta)/u)^2), length(x)-1) < 1-p	
	}
	
	L<-length(x)
	
	if(consistent(x,u,wtm(x,u),p)) {
		return(1:L)
	} else {
		#start search
		found=FALSE
		for( i in 1:(L-1) ) {
			if(verbose) cat("Dropping",i,"\n")

			k<-t(combn(L,L-i))
			
			cset<-apply(k,1,function(ind,x,u,p) consistent(x[ind],u[ind], wtm(x[ind],u[ind]),p),x=x,u=u,p=p)

			if(sum(cset)>0) {				
				found=T
				break
			}
		}
	}
	if(found) {
		k<-k[cset,]
		if(sum(cset)>1) {
			cat(paste("Warning: More than one consistent subset was found in\n     ",format(match.call()),"\n") )		
			if(simplify) {
				FR<-apply(k,1,function(ind,x,u) sum( ((x[ind]-wtm(x[ind],u[ind]))/u[ind] )^2 ), x=x, u=u)
				return( k[which.min(FR),] )
			} 
		}
		return(k)
	} else {
		return(NULL)
	}
}
