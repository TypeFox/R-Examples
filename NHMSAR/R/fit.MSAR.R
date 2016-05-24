fit.MSAR <-
function(
    data,theta,MaxIter=100,eps=1e-5,verbose=FALSE,
covar.emis=NULL,covar.trans=NULL,method=NULL,constraints=FALSE,reduct=FALSE,K=NULL,d.y=NULL,ARfix=FALSE,penalty=FALSE,sigma.diag=FALSE,lambda1=.1,lambda2=.1,a=3.7,...
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
    
    BIC = NULL
    Npar = NULL
    
    cnt <- 0
    FB <- Estep.MSAR(data,theta,covar.trans=covar.trans,covar.emis=covar.emis)
    loglik = FB$loglik
    previous_loglik <- FB$loglik-1000 
    ll_history = NULL
    converged = EM_converged(0,2*eps,eps);
    par = NULL
    while (converged[1]==0 && cnt < MaxIter) {
    	cnt <- cnt+1

    	if (verbose) {print(c("iteration ",cnt,"  loglik = ",loglik),quote = FALSE)}

    	# -----------------------------
    	# ...... M step
    	if (label=='HH') {
    		if (constraints == FALSE & penalty==FALSE & reduct==FALSE) {par = Mstep.hh.MSAR(data,theta,FB) }
    		else if (constraints){
    			par = Mstep.hh.MSAR.with.constraints(data,theta,FB,K=K,d.y=d.y) 
    			attributes(theta)$n_par = M + M*(M-1) + 2*M*d # A and sigma diagonal
    		} 
    		else if (penalty=="ridge")	{
     			par = Mstep.hh.SCAD.cw.MSAR(data,theta,FB,penalty="ridge",lambda1=0,lambda2=lambda2,par=par)
   			#par = Mstep.hh.ridge.MSAR(data,theta,FB,lambda=lambda2)
			}
    		else if (penalty=="LASSO")	{
    			 if (cnt>1) {
    				 par = Mstep.hh.reduct.MSAR(data,theta,FB,sigma.diag=sigma.diag)
    			 } 
    			 else {par = Mstep.hh.lasso.MSAR(data,theta,FB)}
    		}
    		else if (penalty=="SCAD") {
    			par = Mstep.hh.SCAD.MSAR(data,theta,FB,penalty="SCAD",lambda1=lambda1,lambda2=lambda2,par=par)
    		}    
 #   		else if (penalty=="SIS") {
 #   			par = Mstep.hh.SIS.MSAR#(data,theta,FB,penalty="SCAD",lambda1=lambda1,lambda2=lambda2,par=par)
 #   		}
			else if (reduct) {par = Mstep.hh.reduct.MSAR(data,theta,FB,sigma.diag=sigma.diag)}

    		if (order>0) {
    	        theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat)
			} else {
     		    theta=list(par$A0,par$sigma,par$prior,par$transmat)
			} 		}                 
		else if (label=='HN') { 
		    par = Mstep.hn.MSAR(data,theta,FB,covar=covar.emis,verbose=verbose) 
    		if (order>0) {
    			theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat,par$par_emis)
             }       
             else{
             	theta=list(par$A0,par$sigma,par$prior,par$transmat,par$par_emis)
             }
		}
		else if (label=='NH') { 
			par = Mstep.nh.MSAR(data,theta,FB,covar=covar.trans,method=method,ARfix=ARfix,reduct=reduct,sigma.diag=sigma.diag,penalty=penalty,lambda1=lambda1,lambda2=lambda2,par=par) 
    		if (order>0) {
		        theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat,par$par.trans)
             }          
             else{theta=list(par$A0,par$sigma,par$prior,par$transmat,par$par.trans)}
             
             
		}
		else if (label=='NN') { 
			par = Mstep.nn.MSAR(data,theta,FB,covar.emis=covar.emis,covar.trans=covar.trans,method=method) 
    		if (order>0) {
		        theta=list(par$A,par$A0,par$sigma,par$prior,par$transmat,par$par.trans,par$par.emis)
             }
             else{theta=list(par$A0,par$sigma,par$prior,par$transmat,par$par.trans,par$par.emis)}
             
		}
    ll_history[cnt] = loglik
    converged = EM_converged(loglik, previous_loglik, eps)
    previous_loglik = loglik
    attributes(theta)=att.theta
    theta=as.thetaMSAR(theta,label=label,ncov.emis = ncov.emis,ncov.trans=ncov.trans)             

    
        # -----------------------------
    	# ...... E step
    	FB = Estep.MSAR(data,theta,covar.emis=covar.emis,covar.trans=covar.trans)
    	loglik = FB$loglik

 
    }
    #browser()
   	if ( M>1) { 
   		tr = NULL
   		for (m in 1:M) {if (d==1) {tr[m] = theta$sigma[m]}
   			else {tr[m] = sum(diag(theta$sigma[[m]]))}}
   		i.tr = order(tr)
   		theta$A0 = theta$A0[i.tr,]
   		sigma.tmp = NULL
   		A.tmp = theta$A
   		for (m in 1:M) {
   			sigma.tmp[[m]] = theta$sigma[[i.tr[m]]]
   			if (order>0) {
   				if (d>1) {
   					for (o in 1:order ) {
   						A.tmp[[m]][[o]] = theta$A[[i.tr[m]]][[o]]
   					}
   				}
   				else {A.tmp[m,] = theta$A[i.tr[m],]}}

   		}
   		theta$A = A.tmp
   		theta$sigma = sigma.tmp
   		theta$prior = theta$prior[i.tr,]
   		temp = theta$transmat[i.tr,i.tr]
   		attributes(temp)$dimnames[[2]] <- attributes(theta$transmat)$
dimnames[[2]]
		attributes(temp)$dimnames[[2]] <- attributes(theta$transmat)$dimnames[[2]]
		theta$transmat = temp
		if (substr(label,2,2)=="N") {tmp = theta$par.emis 
			for (j in 1:M) {theta$par.emis[[j]] = tmp[[i.tr[j]]]}
		}
		if (substr(label,1,1)=="N") {theta$par.trans = theta$par.trans[i.tr,]} 
		theta = as.thetaMSAR(theta,label=label,ncov.emis = ncov.emis,ncov.trans=ncov.trans)
		FB$probS = FB$probS[,,i.tr]
	}
	#browser()
    if (penalty!="SCAD" ) {lambda1=rep(0,M)}
    npar = M*d+M*(M-1)
    if (substr(label,1,1)=="N") {npar = npar+M*length(theta$par.trans[1,])}
    if (substr(label,2,2)=="N") {npar = npar+M*length(theta$par.emis[1,])}
    for (m in 1:M) {
    	npar = npar+sum(abs(theta$A[[m]][[1]])>0)
    	if (penalty!="SCAD" | max(abs(lambda1))==0) {npar = npar+sum(abs(theta$sigma[[m]][upper.tri(theta$sigma[[m]],diag=TRUE)])>0)}
    	else { npar = npar+sum(abs(par$sigma.inv[[m]][upper.tri(par$sigma.inv[[m]],diag=TRUE)])>1e-5)} # 1e-5 : arbitrary level... 
    }
    attributes(theta)$n_par = npar
	#}
    BIC = -2*ll_history[cnt]+ attributes(theta)$n_par*log(length(c(data)))
    ll.pen = NULL
    if (penalty=="SCAD") {
    	a=3.7
    	if (length(lambda1)==1) {lambda1 = matrix(lambda1,1,M)}
		if (length(lambda2)==1) {lambda2 = matrix(lambda2,1,M)}
    	pen = 0
		for (m in 1:M){
			w = matrix(0,d,d)
			if (lambda1[m]>0) {wi = solve(theta$sigma[[m]])
				abs.S = abs(theta$sigma[[m]])
				w = lambda1[m]*abs.S
				wA = which(abs.S>lambda1[m] & abs.S<=a*lambda1[m])
				w[wA] = -(abs.S[wA]^2-2*a*lambda1[m]*abs.S[wA]+lambda1[m]^2)/2/(a-1)
				wA = which(abs.S>a*lambda1[m])
				w[wA] = (a+1)^2*lambda1[m]^2/2
				w=matrix(w,d,d)
			}
			pen = pen+sum((w-diag(diag(w)))*abs(theta$sigma[[m]]))
			omega = matrix(0,d,d)
			if (lambda2[m]>0) {
				abs.A = abs(theta$A[[m]][[1]])
				omega = lambda2[m]*abs.A
				wA = which(abs.A>lambda2[m] & abs.A<=a*lambda2[m])
				omega[wA] = -(abs.A[wA]^2-2*a*lambda2[m]*abs.A[wA]+lambda2[m]^2)/2/(a-1)
				wA = which(abs.A>a*lambda2[m])
				omega[wA] = (a+1)^2*lambda2[m]^2/2
				omega=matrix(omega,d,d)
			}
			pen = pen+sum(omega*abs(theta$A[[m]][[1]]))
		}
		ll.pen = (FB$loglik-((T-1)*N.samples)*pen)
	}
    res = list(theta=theta,ll_history=ll_history,Iter=cnt,Npar=Npar,BIC=BIC,smoothedprob=FB$probS,ll.pen = ll.pen)
    class(res) <- "MSAR"
    res$call = cl
    res
}
