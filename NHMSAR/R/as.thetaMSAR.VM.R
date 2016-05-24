as.thetaMSAR.VM <-
function(x,label='HH',regime_names=NULL,ncov.emis=ncov.emis,ncov.trans=ncov.trans) {
    if (!is.thetaMSAR(x)) {
        stop('as.thetaMSAR: your input is not like a theta at all')
    }
    
    dimname <- dimnames(x)
    att <- attributes(x)
    att$dimnames <- NULL
    att$names <- NULL
    att$class <- NULL
    
    label = att$label
    M = att$NbRegimes #number of regimes
    order=att$order
    d=att$NbComp
    
    dimname=c("mu","kappa","prior","transmat")
    
    n_par=NULL
    
    if(substr(label,1,1)=="N"){dimname=c(dimname,'par.trans'); n_par=n_par+length(c(x$par.trans))/M}
     if(substr(label,2,2)=="N"){dimname=c(dimname,'par.emis'); n_par=n_par+length(c(x$par.emis))/M}
   
    names(x)=dimname
    
    x$mu = matrix(x$mu,M,d)
    rownames(x$mu)=c(paste("Regime",1:M,sep=""))
    colnames(x$mu)=c(paste("mu",1:d,sep=""))
    
    x$prior = matrix(x$prior,M,d)
    x$prior=matrix(x$prior,M,1)
    rownames(x$prior)=c(paste("Regime",1:M,sep=""))
    colnames(x$prior)=""
    
    rownames(x$transmat)=c(paste("Regime",1:M,sep=""))
    colnames(x$transmat)=c(paste("Regime",1:M,sep=""))
    
    if(substr(label,2,2)=="N"){
    	names(x$par.emis)=c(paste('Regime',1:M,sep=""))
    	for(i in 1:M){
    		if (!is.null(dim(x$par.emis[[i]]))) {rownames(x$par.emis[[i]])=rep("",d)}
     		if (!is.null(dim(x$par.emis[[i]])[2])) {
     			colnames(x$par.emis[[i]])=c(paste('coef.emis',1:dim(x$par.emis[[i]])[2],sep=""))
     		}   	
     	}
    }
    
    if(substr(label,1,1)=="N"){
    	rownames(x$par.trans)=c(paste("Regime",1:M,sep=""))
    	colnames(x$par.trans)=c(paste('coef.trans',1:max(2,ncov.trans+1),sep=""))	
    }



    if( d==1){
    	rownames(x$kappa)=c(paste('Regime',1:M,sep=""))
    	colnames(x$kappa)=c(paste('kappa',1:(order+1),sep=""))
    }
    # else if (order>0){
    	# names(x$A)=c(paste("Regime",1:M,sep=""))
    	# for(i in 1:M){
    		# names(x$A[[i]])=c(paste('A',1:order,sep=""))
        # for(j in 1:order){
          # colnames(x$A[[i]][[j]])=rep("",d)
          # rownames(x$A[[i]][[j]])=rep("",d)
        # }
    		# }
    #}
    
    if (is.null(att$order)) { att$order <- order }
    if (is.null(att$n_par)) { att$n_par <- M^2+M+M*(order+1) }
    
#     if (is.null(att$label)) { att$label <- label } # A VOIR !!!

#    x <- array(x,thetadim,dimnames=dimname)
    for (a in names(att)) {
        attr(x,a) <- att[[a]]
    }
    
    class(x) <- "MSAR.VM"
    return(x)
}
