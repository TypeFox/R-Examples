Mstep.hn.MSAR <-
function(data,theta,FB,covar=NULL,verbose=FALSE)  {  
	N.samples = dim(as.array(data))[2] 
    M = dim(FB$probS)[3]
    order = attributes(theta)$order
	#order = max(p,1)
	T=dim(as.array(data))[1]

    
    d = dim(as.array(data))[3]
    if (is.null(d)|is.na(d)) { d=1 }
    ncov = dim(as.array(covar))[3]
    if (is.null(ncov)|is.na(ncov)) { ncov=1 }
    data = array(data,c(T,N.samples,d)) 
    data2 = array(0,c(order*d,T-order+1,N.samples))
	cpt=1
    for (o in order:1) {
    	for (kd in 1:d) {
    		 data2[cpt,,] = data[o:(T-order+o),,kd]
    		 cpt =cpt+1
    	}
    }
       
	moy = as.matrix(theta$A0,M,d)	
	if(d>1){
	    Sigma = theta$sigma
	    Sigma=matrix(0,M,d^2)
	    for(j in 1:M) { Sigma[j,]=as.vector(theta$sigma[[j]]) }
	    sigma.res = NULL
	}
	else{
		Sigma=theta$sigma
		sigma.res = Sigma
	}
	nh.emissions <- attributes(theta)$nh.emissions 
	#par_emis=matrix(0,M,dim(as.matrix(theta$par.emis[[1]]))[2])
	par_emis=matrix(0,M,length(theta$par.emis[[1]]))
	for (j in 1:M) {
		par_emis[j,] = as.vector(theta$par.emis[[j]]) # Attention a l'ordre dans lequel c'est range...
	}
	
	if(d==1){A <- matrix(0,M,order)}
	else{A <- list()}
	if(d==1){sigma <- matrix(0,M,1)}
	else{sigma <- list()}
	par.emis=list()
	if (order>0) {
	    # A = matrix(0,M,order*d^2)
		# if(d>1){
			# for(j in 1:M){
				# for(i in 1:order){
					# A[j,(((i-1)*d^2+1):(i*d^2))]=as.vector(theta$A[[j]][[i]])
		        # }
		    # }
		# }else{A=theta$A}
		A  = theta$A
		if (attributes(theta)$emis.linear) {
			#T = length(o:(dim(data)[1]-order+o))
			exp_num_trans = 0
			exp_num_visit = 0
			exp_num_visits1 = 0
			postmix = 0
			m = matrix(0,d*order+dim(covar)[3],M) ; m_1 = matrix(0,d*order,M) ; c = matrix(0,M,1) ; s=c ; 
    		op = array(0,c(d*order+dim(covar)[3],d*order+dim(covar)[3],M) ); 
    		op_1 = array(0,c(d*order+dim(covar)[3],d*order,M)) ; 
    		op_2 = array(0,c(d*order,d*order,M) ) ;
			#gamma=matrix(FB$probS[,,j],N.samples,T-order)			
			for (ex in 1:N.samples) { #browser()
				obs = array(data2[,,ex],c(d*order,T))
				xit  = array(FB$probSS[,,,ex],c(M,M,T-2))
				gamma = matrix(FB$probS[ex,1:(T-order),],T-order,M)
				exp_num_trans = exp_num_trans+apply(xit,c(1,2),sum)
				exp_num_visits1 = exp_num_visits1+gamma[1,]
				postmix = postmix+apply(gamma,2,sum)
				obs = t(obs)     					#for (jj in 1:M) {#browser()
				for (j in 1:M) { #browser()
					w = matrix(gamma[,j], 1,T-order)
		  			wobs = cbind(obs[1:(T-order),],covar[1:(T-order),ex,]) * repmat(t(w),1,d*order+dim(covar)[3])
		  			op[,,j] = op[,,j] + t(wobs) %*% cbind(obs[1:(T-order),],covar[1:(T-order),ex,])
          			wobs_1= obs[2:(T-order+1),] * repmat(t(w),1,d*order)
          			if (is.na(sum(obs))) {
            			wobs[is.na(wobs)] = 0
            			wobs_1[is.na(wobs_1)] = 0
            			obs[is.na(obs)] = 0
            		}
            		tmp = apply(wobs,2,sum)
          			m[,j] = m[,j] + tmp
           			m_1[,j] = m_1[,j] + apply(wobs_1, 2,sum)#
            		op_1[,,j] = op_1[,,j] + t(wobs) %*% obs[2:(T-order+1),]
           			op_2[,,j] = op_2[,,j] + t(wobs_1) %*% obs[2:(T-order+1),] #}
    			}	
   			}  
   			for (j in 1:M) { 
			   	Cxx = postmix[j]*op[,,j] - m[,j]%*%t(m[,j])
        	   	Cxy = postmix[j]*op_1[,,j] - m[,j]%*%t(m_1[,j])
        	   	A2tmp = t(Cxy)%*%solve(Cxx) 
        	   	tmp = (m_1[,j]-(A2tmp)%*%m[,j])/postmix[j] 
        	   	moy[j,1:d] = tmp[1:d] 
        	   	#op_1[,,j] = t(as.matrix(op_1[,,j]))
        	   	tmp2 = (op_2[, , j] + (A2tmp) %*% op[, , j] %*% t(A2tmp) - 
                   ((A2tmp) %*% (op_1[, , j]) + t((A2tmp) %*% (op_1[, , j]))))/postmix[j] - tmp %*% t(tmp)
               	if (d==1) {sigma.res[j]=tmp2[1:d,1:d]}
               	else {sigma.res[[j]]=tmp2[1:d,1:d]}
               	dA2 = order*d
               	for(kp in 1:order){
               		#A[j,kp] = A2tmp[1:d,((kp-1)*d+1):(kp*d)] # ok if d==1
               		A[[j]][[kp]] = A2tmp[1:d,((kp-1)*d+1):(kp*d)] # ok if d==1
               		par.emis[[j]] = A2tmp[1:d,(dA2+1):dim(A2tmp)[2]]
               	}
			} 
		} else {
			gamma=matrix(FB$probS[,,j],N.samples,T-order)
			for (j in 1:M){	    A2 = matrix(0,M,order*d^2)
		 	if(d>1){
			 	for(j in 1:M){
					for(i in 1:order){
						A2[j,(((i-1)*d^2+1):(i*d^2))]=as.vector(theta$A[[j]][[i]])
		         	}
		     	}
		 	}else{ A2=theta$A }
			if (d==1) {a = as.matrix(c(A2[j,],moy[j,],log(Sigma[j,]),par_emis[j,]))}
			else {a = as.matrix(c(A2[j,],moy[j,],Sigma[j,],par_emis[j,]))}
            resopt = optim(a,ll_gausshn.MSAR,data=data,gamma=gamma,order=order,nh.emissions=nh.emissions,covar=covar, gr = NULL, method = "L-BFGS-B",control = list(trace=TRUE), hessian = FALSE)      
         	
            if(d==1){ A[j,]=resopt$par[1:order]
        	}else{ A[[j]]=list()
            	for(o in 1:order){ A[[j]][[o]]=matrix(resopt$par[((o-1)*d^2+1):o*d^2],d,d)}
        	}
        	moy[j,]=resopt$par[(order*d^2+1):(order*d^2+d)]
        	if(d==1){ sigma.res[j,] <-exp(resopt$par[order+2]) #Sigma[j,] <-(resopt$par[order+2])
        	}else { sigma.res[[j]]=matrix(resopt$par[(order*d^2+d+1):((order+1)*d^2+d)],d,d) }
			par.emis[[j]]=matrix(resopt$par[((order+1)*d^2+d+1):length(resopt$par)],d,ncov)
		}
      }	
    }	
     else { # p (order) = 0
    	for (j in 1:M) {
    		a = as.matrix(c(moy[j,],log(Sigma[j,]),par_emis[j,]))
   		    gamma=matrix(FB$probS[,,j],N.samples,T)
    		resopt = ucminf(a,fn=ll_gausshn_p0.MSAR,data=data[(order+1):T,,],gamma=gamma,nh.emissions=nh.emissions,covar=array(covar[(order+1):T,,],c(T-order,N.samples,dim(covar)[3])), gr = NULL, control = list(trace=FALSE))          	
            moy[j,]=resopt$par[1:d]
            if(d==1){sigma.res[j,] <- exp(resopt$par[2])}               
            else{sigma.res[[j]] <- matrix(exp(resopt$par[(d+1):(d+d^2)]))}
            par.emis[[j]]=matrix(resopt$par[(d+d^2+1):length(resopt$par)],d,ncov)
        }
	}
	exp_num_trans = 0
	exp_num_visits1 = 0    
	for (ex in 1:N.samples) {
		xit  = array(FB$probSS[,,,ex],c(M,M,T-2)) # xit = FB$probSS[,,,ex]
		gamma = matrix(FB$probS[ex,1:(T-order),],T-order,M) # gamma = FB$probS[ex,,]
		exp_num_trans = exp_num_trans+apply(xit,c(1,2),sum) # sum3(xit) = sum_{t=1}^T P(S_t=i,S_{t+1}=j|y_1,...,y_T
		exp_num_visits1 = exp_num_visits1+gamma[1,] # gamma[,1] = P(S_1=i|y_1,...,y_T)		    

	}	        
    prior = normalise(exp_num_visits1)
    transmat = mk_stochastic(exp_num_trans)
    if (order>0) {
	    list(A=A,sigma=sigma.res,A0=moy,prior=prior,transmat=transmat,par_emis=par.emis)
	} else {
		tmp = list(sigma=sigma.res,A0=moy,prior=prior,transmat=transmat,par_emis=par.emis)
		return(tmp)
	}
}
