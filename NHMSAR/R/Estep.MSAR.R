Estep.MSAR <-
function(
    data,theta,smth=FALSE,verbose=FALSE,
    covar.emis=covar.emis,covar.trans=covar.trans
) {
	go <- Sys.time()
    if (verbose) { message('Starting LWS_pfilt') }
    att.theta = attributes(theta)
    label <- att.theta$label
    p <- att.theta$order
    if (verbose) print(theta)
    M <- att.theta$NbRegimes
    d <- att.theta$NbComp
    if(is.null(d) || is.na(d)){d=1}
    n_par <- att.theta$n_par
    order = att.theta$order
   
    data <- as.array(data)
    T = dim(data)[1]
    if (is.null(T)) {T = length(data)}
    if (abs(length(data)/T-trunc(length(data)/T)) >1e-5) {stop('error : size of data should be nxT.sample with n integer')}
    N.samples = dim(data)[2]
    if(is.null(N.samples) || is.na(N.samples)){N.samples <- 1}
    data <- array(data,c(T,N.samples,d))
    
    gamma <- array(0,c(N.samples,T-p,M))
    xi = array(0,c(M,M,T-(p+1),N.samples))
    loglik = 0
    if(substr(label,2,2) == 'N'){
    	ncov.emis = dim(covar.emis)[3]
    	if(is.null(ncov.emis) || is.na(ncov.emis)){ncov.emis=1}
    	covar.emis=array(covar.emis,c(T,N.samples,ncov.emis))
    }
    if (substr(label,1,1) == 'H') { 
    	tmat = as.matrix(theta$transmat)
		prior = as.matrix(theta$prior)
    	for (ex in 1:N.samples) {
    		if (verbose) { print(ex) }
   			g <- emisprob.MSAR(data[,ex,],theta=theta,covar=covar.emis[,ex,])
   			FB = forwards_backwards(prior, tmat, g)
   			gamma[ex,,] = t(FB$gamma)
   			xi[,,,ex] = FB$xi
   			loglik = loglik + FB$loglik
    	   	}
     }
     else  {
        if(missing(covar.trans)){stop("error : covariable is missing")}
        if (length(covar.trans)==1) {
    		Lag = covar.trans+1
    		covar.trans = array(data[(1):(T-Lag+1),,],c(T-Lag+1,N.samples,d))
    		data =  array(data[Lag:T,,],c(T-Lag+1,N.samples,d))
    		#T = dim(data)[1]
    	}
		ncov.trans = dim(covar.trans)[3]
    	if(is.null(ncov.trans) || is.na(ncov.trans)){ncov.trans=1}
    	ct = array(0,c(T,N.samples,ncov.trans))
    	ct[1:min(dim(ct)[1],dim(covar.trans)[1]),,] = covar.trans[1:min(dim(ct)[1],dim(covar.trans)[1]),,]
		covar.trans=ct
     	for (ex in 1:N.samples) {
    		if (verbose) { print(ex) }
    		g <- emisprob.MSAR(data[,ex,],theta=theta,covar=covar.emis[,ex,])
    		transmat = theta$transmat
    		par.trans = theta$par.trans
    		nh_transition = attributes(theta)$nh.transitions
     		inp = covar.trans[(order+1):T,ex,]
    		transmat.t = nh_transition(array(inp,c((T-order),1,ncov.trans)),par.trans,transmat)
            FB = nhforwards_backwards(theta$prior, transmat.t, g)
     		gamma[ex,,] = t(FB$gamma)
    	    xi[,,,ex] = FB$xi
    	    loglik = loglik + FB$loglik
         } 
    }
    list(loglik=loglik,probS=gamma,probSS=xi,M=FB$M)
}
