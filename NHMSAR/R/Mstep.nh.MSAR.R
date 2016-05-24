Mstep.nh.MSAR <-
function(data,theta,FB,covar=NULL,method=method,ARfix=FALSE,reduct=FALSE,penalty=FALSE,sigma.diag=FALSE,lambda1=lambda1,lambda2=lambda2,par = NULL)  {  
	
	  order=attributes(theta)$order
	  d=dim(data)[3]
	  if(is.na(d) | is.null(d)){d=1}
	  M=attributes(theta)$NbRegimes
	  
      if (length(covar)==1) {
    	  Lag = covar
    	  covar = array(data[(1):(T-Lag+1),,],c(T-Lag+1,N.samples,d))
    	  data =  array(data[Lag:T,,],c(T-Lag+1,N.samples,d))
      }

	  N.samples = dim(covar)[2]
	  ncov.trans = dim(covar)[3]
      if (!ARfix & !reduct & penalty==FALSE) {par.hh = Mstep.hh.MSAR(data,theta,FB)}
      else if (reduct) {par.hh = Mstep.hh.reduct.MSAR(data,theta,FB,sigma.diag=sigma.diag)}
      else if (penalty=="SCAD") { par.hh = Mstep.hh.SCAD.MSAR(data,theta,FB,penalty=penalty,lambda1=lambda1,lambda2=lambda2,par=par)
    		}
      else {
      	par.hh = NULL
      	par.hh$A=theta$A
      	par.hh$A0=theta$A0
      	par.hh$sigma=theta$sigma
      	par.hh$prior=theta$prior
      }
      theta$transmat[which(theta$transmat<1e-15)] = 1e-15
      theta$transmat = mk_stochastic(theta$transmat)
      trans = para_trans(theta$transmat) 
      par.trans = theta$par.trans
      nh_transition = attributes(theta)$nh.transitions
      par.init = plie2(trans,par.trans)   
 	  lxi = dim(FB$probSS)[3]

      if (order>0) {deb = order+1}
      else {deb = 1} 
      if (is.null(method)) {method="ucminf"}
      if (method=="ucminf") {resopt = ucminf(par.init,fn=loglik_nh_inp,gr=NULL,covar=array(covar[deb+(1:(lxi)),,],c(lxi,N.samples,ncov.trans)),xi=FB$probSS,nh_transition=nh_transition,hessian=0,control = list(trace=FALSE))}
      else if (method=="L-BFGS-B"){resopt = optim(par.init,fn=loglik_nh_inp,gr=NULL,covar=array(covar[deb+(1:(lxi)),,],c(lxi,N.samples,ncov.trans)),xi=FB$probSS,nh_transition=nh_transition,hessian=0,control = list(trace=FALSE),method="L-BFGS-B")}
	  else if (method=="BFGS"){resopt = optim(par.init,fn=loglik_nh_inp,gr=NULL,covar=array(covar[deb+(1:(lxi)),,],c(lxi,N.samples,ncov.trans)),xi=FB$probSS,nh_transition=nh_transition,hessian=0,control = list(trace=FALSE),method="BFGS")}
      res = deplie2(resopt$par);
      trans = res$trans ;
      par.trans = res$par ;
      transmat=para_trans_inv(trans);

     if (order>0) {
	    list(A=par.hh$A,sigma=par.hh$sigma,A0=par.hh$A0,prior=par.hh$prior,transmat=transmat,par.trans=par.trans,sigma.inv=par.hh$sigma.inv)
	 } else {
		list(sigma=par.hh$sigma,A0=par.hh$A0,prior=par.hh$prior,transmat=transmat,par.trans=par.trans,sigma.inv=par.hh$sigma.inv)
	 }
}
