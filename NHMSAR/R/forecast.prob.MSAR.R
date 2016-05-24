forecast.prob.MSAR <-
function(data,theta,yrange=NULL,covar.emis=NULL,covar.trans=NULL){
	# yrange : only 1d???
	if (is.null(yrange)) {
		yrange = seq(min(data,na.rm=TRUE),max(data,na.rm=TRUE),length.out=100)
	}
	lyr = length(yrange)
	T = dim(data)[1]
	N.samples = dim(data)[2]
	d = dim(data)[3]
	M = attributes(theta)$NbRegimes
	label = attributes(theta)$label
	order <- attributes(theta)$order
	if (substr(label,1,1) == "N") {
		if (length(covar.trans)==1) {
    		Lag = covar.trans+1
    		covar.trans = array(data[(1):(T-Lag+1),,],c(T-Lag+1,N.samples,d))
    		data =  array(data[Lag:T,,],c(T-Lag+1,N.samples,d))
    	}
    	ncov.trans = dim(covar.trans)[3]
    	if(is.null(ncov.trans) || is.na(ncov.trans)){ncov.trans=1}
		covar.trans=array(covar.trans,c(T,N.samples,ncov.trans))
	}
	prob = array(0,c(lyr,T,N.samples))
	Yhat = matrix(0,T,N.samples)
	Yhat0 = Yhat
	for (ex in 1:N.samples){
		for (o in (1:order)) {
			for (m in 1:M) {
				prob[,o,ex] = prob[,o,ex]+pdf.norm(matrix(yrange,1,lyr),theta$A0[m],matrix(theta$sigma[[m]]))*theta$prior[m]
			}
		   	y = yrange*prob[,o,ex]
		   	Yhat[o,ex] = sum(diff(yrange)*(y[-length(y)]+y[-1])/2) # sam as trapz of matlab
		}
		gy.tmp <- emis4ps.MSAR(yrange,data[,ex,],theta,covar=covar.emis[,ex,])
		gy = array(0,c(lyr,T,M))
		gy[,(order+1):T,] = gy.tmp
   		g = emisprob.MSAR(data[,ex,],theta,covar=covar.emis[,ex,])
		prior = as.matrix(theta$prior)
		tmat = as.matrix(theta$transmat)
		if (substr(label,1,1) == "H") {
  			FB = forwards_backwards(prior, tmat, g) 
  			transmat.t = array(tmat,c(M,M,T))
  		} else {
    		par.trans = theta$par.trans
    		nh_transition = attributes(theta)$nh.transitions
    		inp = covar.trans[(order+1):T,ex,]
    		transmat.t = nh_transition(array(inp,c((T-order),1,ncov.trans)),par.trans,tmat)
            FB = nhforwards_backwards(theta$prior, transmat.t, g)
            trans.tmp = array(0,c(M,M,T))
            trans.tmp[,,(order+1):T] = transmat.t
            transmat.t = trans.tmp
        }

   		alpha = matrix(0,M,T)
   		alpha[,1] = prior
   		alpha[,(order+1):T] = FB$alpha		
   		w = matrix(NA,M,T)
   		for (t in (order+2):T) {
   			w[,t] = (transmat.t[,,t]%*%alpha[,t-1]) # alpha[,1] = prior * obslik[,t] et t=order+1
   			w[,t] = w[,t]/sum(w[,t])
   			prob[,t,ex] = gy[,t,]%*%w[,t]
   			y = yrange*prob[,t,ex]
   			Yhat[t,ex] = sum(diff(yrange)*(y[-length(y)]+y[-1])/2) # sum as trapz of matlab
   			y = yrange*(gy[,t,]%*%alpha[,t])
   			Yhat0[t,ex] = sum(diff(yrange)*(y[-length(y)]+y[-1])/2) # sum as trapz of matlab
   		}	
   	}
	list(yrange=yrange,prob=prob,Yhat=Yhat,Yhat0=Yhat0)
}
