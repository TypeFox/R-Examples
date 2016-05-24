fit.MSAR.VM <-
function(
    data,theta,MaxIter=100,eps=1e-5,verbose=FALSE,
    covar.emis=NULL,covar.trans=NULL,method=NULL,constr=0,...
) { 
	cl <- match.call()
    now <- Sys.time()
    if (missing(theta)) {
    	stop('can not fit a MSAR model without initial value for theta')
    }
    
    att <- attributes(data)   
    att.theta <- attributes(theta)
    
    data <- as.array(data)
    
    T <- dim(data)[1]
    if(is.null(T) || is.na(T)){T <- length(data)}
    N.samples <- dim(data)[2]
    if(is.null(N.samples) || is.na(N.samples)){N.samples <- 1}
    d <- att.theta$NbComp
    if(is.null(d)|is.na(d)){d <- 1}
    
    data <- array(data,c(T,N.samples,d))
    
    order <- att.theta$order
    label <- att.theta$label
    M <- att.theta$NbRegimes
    
    if ((missing(covar.trans) && substr(label,1,1)=="N") || (missing(covar.emis) && substr(label,2,2)=="N")) {
    	print("Can not fit a non homogeneous MSAR without covariable")
    }
    if (length(covar.trans)==1) {
    	Lag = covar.trans+1
    	covar.trans = array(data[(1):(T-Lag+1),,],c(T-Lag+1,N.samples,d))
    	data =  array(data[Lag:T,,],c(T-Lag+1,N.samples,d))
    }

    if (!is.null(covar.emis)) {
    	if (is.null(dim(covar.emis)) ) {ncov.emis=1}
    	else if (is.na(dim(covar.emis)[3])) {ncov.emis=1}
    	else if (!is.na(dim(covar.emis)[3])) {ncov.emis = dim(covar.emis)[3]}
    	covar.emis = array(covar.emis,c(T,N.samples,ncov.emis))
    }
    		# si dim(covar.emis)[3] = NA, mattre sous forme d'array
    else {ncov.emis = 0}
	if (!is.null(covar.trans)) {ncov.trans = dim(covar.trans)[3]}
    else {ncov.trans = 0}
    #if (substr(label,1,1)=="N") {covar = as.matrix(covar,T,length(covar)/T)}
    
    BIC = NULL
    Npar = NULL
    
    cnt <- 0
    FB <- Estep.MSAR.VM(data,theta,covar.trans=covar.trans,covar.emis=covar.emis)
    loglik = FB$loglik
    previous_loglik <- FB$loglik-1000 
    ll_history = NULL
    converged = EM_converged(0,2*eps,eps);
    while (converged[1]==0 && cnt < MaxIter) {
    	cnt <- cnt+1

    	if (verbose) {print(c("iteration ",cnt,"  loglik = ",loglik),quote = FALSE)}
 
    	# -----------------------------
    	# ...... M step
    	if (label=='HH') {
    		par = Mstep.hh.MSAR.VM(data,theta,FB,constr=constr)  		
     	theta=list(par$mu,par$kappa,par$prior,par$transmat)
		Npar = M*(order+2)+M-1+M*(M-1)
	} # else if (label=='HN') { 
		# par = Mstep.hn.MSAR(data,theta,FB,covar=covar.emis,verbose=verbose) 
    		# if (order>0) {
    			# theta=list(par$A,par$mu,par$sigma,par$prior,par$transmat,par$par_emis)
         # }else{
             	# theta=list(par$mu,par$sigma,par$prior,par$transmat,par$par_emis)
             # }
             # Npar = M*(order+2)+M-1+M*(M-1)+length(c(par$par.emis))                                         
		# }
	else if (label=='NH') { 
		par = Mstep.nh.MSAR.VM(data,theta,FB,covar.trans=covar.trans,method=method,constr=constr) 
    		theta = list(par$mu,par$kappa,par$prior,par$transmat,par$par.trans)
         Npar = M*(order+2)+M-1+M*(M-1)+M
	}
		# else if (label=='NN') { 
			# par = Mstep.nn.MSAR(data,theta,FB,covar.emis=covar.emis,covar.trans=covar.trans,method=method) 
    		# if (order>0) {
		        # theta=list(par$A,par$mu,par$sigma,par$prior,par$transmat,par$par.trans,par$par.emis)
             # }
             # else{theta=list(par$mu,par$sigma,par$prior,par$transmat,par$par.trans,par$par.emis)}
             
             # Npar = M*(order+2)+M-1+M*(M-1)+length(c(par$par.emis))+length(c(par$par.trans))
		# }
    ll_history[cnt] = loglik
    converged = EM_converged(loglik, previous_loglik, eps)
    previous_loglik = loglik
    attributes(theta)=att.theta
    theta=as.thetaMSAR.VM(theta,label=label,ncov.emis = ncov.emis,ncov.trans=ncov.trans)       

    
        # -----------------------------
    	# ...... E step
    	FB = Estep.MSAR.VM(data,theta,covar.emis=covar.emis,covar.trans=covar.trans)
    	loglik = FB$loglik


    }
   BIC = -2*ll_history[cnt]+ Npar*log(length(c(data)))
    res = list(theta=theta,ll_history=ll_history,Iter=cnt,Npar=Npar,BIC=BIC,smoothedprob=FB$probS)
    class(res) <- "MSAR"
    res$call = cl
    res
}
