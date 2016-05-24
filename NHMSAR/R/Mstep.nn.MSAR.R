Mstep.nn.MSAR <-
function(data,theta,FB,covar.trans=covar.trans,covar.emis=covar.emis,method=NULL)  {  
	
	  order=attributes(theta)$order
	  d=dim(data)[3]
	  if(is.na(d) | is.null(d)){d=1}
	  M=attributes(theta)$NbRegimes
	  
      if (length(covar.trans)==1) {
    	  Lag = covar.trans
    	  covar.trans = array(data[(1):(T-Lag+1),,],c(T-Lag+1,N.samples,d))
    	  data =  array(data[Lag:T,,],c(T-Lag+1,N.samples,d))
      }

      T = dim(data)[1]
      N.samples = dim(data)[2]
      ncov.trans = dim(covar.trans)[3]
     
	  # Lcov <- dim(covar.trans)[1]
	  # if(is.null(Lcov) || is.na(Lcov)){Lcov <- length(covar.trans)}
	  # lcov <- dim(covar)[2]
	  # if(is.null(lcov) || is.na(lcov)){lcov <- 1}
	  #covar <- matrix(covar,Lcov,lcov)
	  
      par.hn = Mstep.hn.MSAR(data,theta,FB,covar.emis)
      theta$transmat[which(theta$transmat<1e-15)] = 1e-15
      theta$transmat = mk_stochastic(theta$transmat)
      trans = para_trans(theta$transmat) 
      par.trans = theta$par.trans
      nh_transition = attributes(theta)$nh.transitions
      par.init = plie2(trans,par.trans)   
      lxi = dim(FB$probSS)[3]
       if (order>0) {deb = order+1}
      else {deb = 1} 
 resopt = ucminf(par.init,fn=loglik_nh_inp,gr=NULL,covar=array(covar.trans[deb+(1:(lxi)),,],c(lxi,N.samples,ncov.trans)),xi=FB$probSS,nh_transition=nh_transition,hessian=0,control = list(trace=FALSE,xtol=1e-4))
      res = deplie2(resopt$par,M);
      trans = res$trans ;
      par.trans = res$par ;
      transmat=para_trans_inv(trans);

     if (order>0) {
	    list(A=par.hn$A,sigma=par.hn$sigma,A0=par.hn$A0,prior=par.hn$prior,par.emis=par.hn$par_emis,transmat=transmat,par.trans=par.trans)
	 } else {
		list(sigma=par.hn$sigma,A0=par.hn$A0,prior=par.hn$prior,par.emis=par.hn$par_emis,transmat=transmat,par.trans=par.trans)
	 }
}
