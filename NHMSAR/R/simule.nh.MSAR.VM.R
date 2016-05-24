simule.nh.MSAR.VM <- function(theta,Y0,T,N.samples = 1,covar.emis=NULL,covar.trans=NULL) {
	
	M = attributes(theta)$NbRegimes
	order = attributes(theta)$order
	d <- attributes(theta)$NbComp
	if (length(Y0) < order) {stop(paste("Length of Y0 should be equal to ",order,sep=""))}
	label = attributes(theta)$label
	
	L = 0
	Y = array(0,c(T,N.samples,d))
	S = matrix(0,T,N.samples)
	Y[1:max(order,1),,] = Y0[1:max(order,1),,]
	transition =  array(0,c(M,M,T,N.samples))
	
	for (ex in 1:N.samples) {
		S[1,ex] = which(rmultinom(1, size = 1, prob = theta$prior)==1)
	}
	
	if (substr(label,1,1) == "N") {
		nh_transitions = attributes(theta)$nh.transitions
		if (length(covar.trans)==1) {
			L = covar.trans
		} else {
#transition=nh_transitions(c(covar.trans),theta$par.trans,theta$transmat);
			for (ex in 1:N.samples) {transition[,,,ex]=nh_transitions(array(covar.trans[,ex,],c(T,1,dim(covar.trans)[3])),theta$par.trans,theta$transmat);}
		}
	} else {
		transition = array(theta$transmat,c(M,M,T,N.samples))
	}
	
	if (order>1) {
		for (ex in 1:N.samples) {
			for (t in 2:order){ 
				S[t,ex] = which.max(rmultinom(1, size = 1, prob = theta$transmat[S[t-1,ex],]))
			}
		}
	}
	kappa = theta$kappa
	for (ex in 1:N.samples) {
		if (order>0) {
			for(t in max(c(2,order+1)):(T)){
				if (substr(label,1,1) == "N" & length(covar.trans)==1) {
					transition[,,t,ex] = nh_transitions(array(Y[t-L,ex,],c(1,1,dim(Y0)[3])),theta$par.trans,theta$transmat)
				}
				S[t,ex] = which.max(rmultinom(1, size = 1, prob = transition[S[t-1,ex], ,t,ex]))
				
				para_com = kappa[S[t,ex],1]*exp(1i*theta$mu[S[t,ex]])
				for (o in 1:order) {
					para_com = para_com+kappa[S[t,ex],o+1]*exp(1i*Y[ t-order+o-1,ex,1])
				}
#sum(kappa[S[t,ex],2:length(kappa[1,])]*exp(1i*Y[(t-order):(t-1),ex,1]))
				Y[t,ex,1] = sample.VM(Arg(para_com),Mod(para_com))
			}
		} else {
			for (t in max(c(2,order+1)):T) {
				if (substr(label,1,1) == "N" & length(covar.trans)==1) {
					transition[,,t,ex] = nh_transitions(matrix(Y[t-L,ex,],1,1),theta$par.trans,theta$transmat)
				}
				S[t,ex] = which.max(rmultinom(1, size = 1, prob = transition[S[t-1,ex], ,t,ex]))
				
				para_com = kappa[S[t,ex],1]*exp(1i*theta$mu[S[t,ex]])
				Y[t,ex,1] = sample.VM(Arg(para_com),Mod(para_com))
			}
		}
	}
	return(list(S=S,Y=Y))
}

