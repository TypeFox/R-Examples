Mstep.hh.MSAR.with.constraints <-
function(data,theta,FB,K,d.y)  {  
	T=dim(data)[1]
	N.samples = dim(as.array(data))[2] 
	d = dim(as.array(data))[3] 
	if(is.null(d)|is.na(d)) {d = 1}
    M <- attributes(theta)$NbRegimes
    p <- attributes(theta)$order
	order <- max(p,1)
	data = array(data,c(T,N.samples,d))
	
	# Markov chain
	T.tmp = length(1:(T-order+1))
	exp_num_trans = 0
	exp_num_visit = 0
	exp_num_visits1 = 0
	for (ex in 1:N.samples) {
		xit  = array(FB$probSS[,,,ex],c(M,M,T.tmp-2))
		gamma = matrix(FB$probS[ex,1:(T.tmp-1),],T.tmp-1,M)
		exp_num_trans = exp_num_trans+apply(xit,c(1,2),sum)
		exp_num_visits1 = exp_num_visits1+gamma[1,]
	}			       
    prior = normalise(exp_num_visits1)
    prior=matrix(prior,M,1)
    transmat = mk_stochastic(exp_num_trans)
    
    # AR part
    A2.tot = list()
    sigma = list()
    moy = matrix(0,M,d)
    for (j in 1:M) {
    	A2.tot[[j]] = list()
    	for (kp in 1:order) {A2.tot[[j]][[kp]] = matrix(0,d,d)}
    	sigma[[j]] = matrix(0,d,d)
    }
    A2.tot[[j]][[kp]] = matrix(0,d,d)
    ik = NULL
	for (k in 1:K) {
		for ( i in 1:d.y) {ik[i] = K*(i-1)+k }
		data.tmp = array(data[,,ik],c(T,N.samples,d.y))
		#data.tmp = data[,,ik]
		data2 = array(0,c(order*d.y,T-order+1,N.samples))
		cpt=1
    	for (o in order:1) {
    		for (kd in 1:d.y) {
    		 	data2[cpt,,] = data.tmp[o:(T-order+o),,kd]
    		 	cpt =cpt+1
    		}
    	}
        
		postmix = 0
		m = matrix(0,d.y*order,M) ; m_1 = m ; c = matrix(0,M,1) ; s=c ; 
    	op = array(0,c(d.y*order,d.y*order,M) ); op_1 = op ; op_2 = op_1 ;

		for (ex in 1:N.samples) {
			obs = array(data2[,,ex],c(d.y*order,T.tmp))
			xit  = array(FB$probSS[,,,ex],c(M,M,T.tmp-2))
			gamma = matrix(FB$probS[ex,1:(T.tmp-1),],T.tmp-1,M)
			postmix = postmix+apply(gamma,2,sum)
			obs = t(obs)     
			for (j in 1:M) {
		   		w = matrix(gamma[,j], 1,T.tmp-1)
		    	wobs = obs[1:(T.tmp-1),] * repmat(t(w),1,d.y*order)
            	wobs_1= obs[2:T.tmp,] * repmat(t(w),1,d.y*order)
            	if (is.na(sum(obs))) {
            		wobs[is.na(wobs)] = 0
            		wobs_1[is.na(wobs_1)] = 0
            		obs[is.na(obs)] = 0
            	}            	
            	m[,j] = m[,j] + apply(wobs,2,sum)
            	m_1[,j] = m_1[,j] + apply(wobs_1, 2,sum)
            	op[,,j] = op[,,j] + t(wobs) %*% obs[1:(T.tmp-1),]
            	op_1[,,j] = op_1[,,j] + t(wobs) %*% obs[2:T.tmp,]
            	op_2[,,j] = op_2[,,j] + t(wobs_1) %*% obs[2:T.tmp,]
         	}
     	}      
        if (min(postmix)<1e-6) {stop("error : smoothing probabilities are to small, in one regime at least. You should revise initialisation.")}
        #moy <- array(0,c(M,d.y))
        Sigma <- array(0,c(d.y,d.y,M))
		A2 <-array(0,c(M,order,d.y^2))
		if (p>0) {																		
		   for (j in 1:M) { 
			   Cxx = postmix[j]*op[,,j] - m[,j]%*%t(m[,j])
        	   Cxy = postmix[j]*op_1[,,j] - m[,j]%*%t(m_1[,j])
        	   A2tmp = t(Cxy)%*%solve(Cxx) 
        	   tmp = (m_1[,j]-(A2tmp)%*%m[,j])/postmix[j] 
        	   moy[j,ik] = tmp[1:d.y] 
        	   op_1[,,j] = t(as.matrix(op_1[,,j]))
        	   tmp2 = (op_2[, , j] + (A2tmp) %*% op[, , j] %*% t(A2tmp) - 
                   ((A2tmp) %*% t(op_1[, , j]) + t((A2tmp) %*% t(op_1[, , j]))))/postmix[j] - tmp %*% t(tmp)       
               for (ii in 1:d.y) {
					for (jj in 1:d.y) {
						for (kp in 1:order){ 
							A2.tot[[j]][[kp]][ik[ii],ik[jj]] =  matrix(A2tmp[1:d.y,((kp-1)*d.y+1):(kp*d.y)],d.y,d.y)[ii,jj]
					}
					sigma[[j]][ik[ii],ik[jj]] = tmp2[ii,jj]
				}
				}
			}

		} else {
			for (j in 1:M) { 
		    	tmp = (m[,j])/postmix[j] 
			    moy[j,ik] = tmp[1:d.y] 
			    tmp2 = op[, , j]/postmix[j] - tmp %*% t(tmp)
                for (i in 1:d.y) {
                	for (j in 1:d.y) {
                		sigma[[j]][ik[ii],ik[jj]] = tmp2[ii,jj]
                	}
                }
			}
		}
	}	
 																
    if (p>0) {
	    list(A=A2.tot,A0=moy,sigma=sigma,prior=prior,transmat=transmat)
	} else {
		list(sigma=sigma,A0=moy,prior=prior,transmat=transmat)
	}
}
