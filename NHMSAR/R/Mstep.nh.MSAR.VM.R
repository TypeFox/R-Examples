Mstep.nh.MSAR.VM <-
function(data,theta,FB,covar.trans=NULL,method=method,constr=0)  {  
	
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
      par.hh = Mstep.hh.MSAR.VM(data,theta,FB,constr)
      theta$transmat[which(theta$transmat<1e-15)] = 1e-15
      theta$transmat = mk_stochastic(theta$transmat)
      trans = para_trans(theta$transmat) 
      par.trans = theta$par.trans
      nh_transition = attributes(theta)$nh.transitions
      #par.init = plie2(trans,par.trans)   
      par.init = plie2.VM(trans,par.trans)   
 	  lxi = dim(FB$probSS)[3]
 
      if (order>0) {deb = order+1}
      else {deb = 1} 
	resopt = ucminf(par.init,fn=loglik_nh_inp.VM,gr=NULL,covar=array(covar[deb+(1:(lxi)),,],c(lxi,N.samples,ncov.trans)),xi=FB$probSS,nh_transition=nh_transition,hessian=0,control = list(trace=FALSE))
      #res = deplie2(resopt$par);
      res = deplie2.VM(resopt$par);
      trans = res$trans ;
      par.trans = res$par ;
      transmat=para_trans_inv(trans);

	    list(mu=par.hh$mu,kappa = par.hh$kappa,prior=par.hh$prior,transmat=transmat,par.trans=par.trans)
}
