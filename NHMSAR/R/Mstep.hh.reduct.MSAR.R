Mstep.hh.reduct.MSAR <-
function(data,theta,FB,sigma.diag=FALSE)  {  
	T=dim(data)[1]
	N.samples = dim(as.array(data))[2] 
	d = dim(as.array(data))[3]
	if(is.null(d)|is.na(d)) {d = 1}
    M <- attributes(theta)$NbRegimes
    p <- attributes(theta)$order
	order <- max(p,1)
	data = array(data,c(T,N.samples,d))
	data2 = array(0,c(order*d,T-order+1,N.samples))
	cpt=1
    for (o in order:1) {
    	for (kd in 1:d) {
    		 data2[cpt,,] = data[o:(T-order+o),,kd]
    		 cpt =cpt+1
    	}
    }
    T = length(o:(dim(data)[1]-order+o))
	exp_num_trans = 0
	exp_num_visit = 0
	exp_num_visits1 = 0
	postmix = 0
	m = matrix(0,d*order,M) ; m_1 = m ; c = matrix(0,M,1) ; s=c ; 
    op = array(0,c(d*order,d*order,M) ); op_1 = op ; op_2 = op_1 ;
	for (ex in 1:N.samples) {
		obs = array(data2[,,ex],c(d*order,T))
		xit  = array(FB$probSS[,,,ex],c(M,M,T-2))
		gamma = matrix(FB$probS[ex,1:(T-1),],T-1,M)
		exp_num_trans = exp_num_trans+apply(xit,c(1,2),sum)
		exp_num_visits1 = exp_num_visits1+gamma[1,]
		postmix = postmix+apply(gamma,2,sum)
		obs = t(obs)   
		
		for (j in 1:M) {
		   	w = matrix(gamma[,j], 1,T-1)
 		    wobs = obs[1:(T-1),] * repmat(t(w),1,d*order)
            wobs_1= obs[2:T,] * repmat(t(w),1,d*order)
            if (is.na(sum(obs))) {
            	wobs[is.na(wobs)] = 0
            	wobs_1[is.na(wobs_1)] = 0
            	obs[is.na(obs)] = 0
            }
            m[,j] = m[,j] + apply(wobs,2,sum)
            m_1[,j] = m_1[,j] + apply(wobs_1, 2,sum)
            op[,,j] = op[,,j] + t(wobs) %*% obs[1:(T-1),]
            op_1[,,j] = op_1[,,j] + t(wobs) %*% obs[2:T,]
            op_2[,,j] = op_2[,,j] + t(wobs_1) %*% obs[2:T,]
         }
    }

    # Markov chain
    prior = normalise(exp_num_visits1)
    prior=matrix(prior,M,1)
    transmat = mk_stochastic(exp_num_trans)
        
    if (min(postmix)<1e-6) {stop("error : smoothing probabilities are to small, in one regime at least. You should revise initialisation.")}
    moy <- array(0,c(M,d))
    sigma <- list()
	A2 <-list()
																		
	for (j in 1:M) {
		#S = theta$sigma[[j]]
		#Si = solve(S)
		A2[[j]] = list()
		# Cxx = postmix[j]*op[,,j] - m[,j]%*%t(m[,j])
        # Cxy = postmix[j]*op_1[,,j] - m[,j]%*%t(m_1[,j])
  		Cxx = (postmix[j]*op[,,j] - m[,j]%*%t(m[,j]))/postmix[j]^2
        Cxy = (postmix[j]*op_1[,,j] - m[,j]%*%t(m_1[,j]))/postmix[j]^2
        Cyy = (postmix[j]*op_2[,,j] - m_1[,j]%*%t(m_1[,j]))/postmix[j]^2

        A.th = theta$A[[j]][[1]]
        S.th = theta$sigma[[j]]
		SA.th = Cyy-Cxy%*%A.th-A.th%*%Cxy+A.th%*%Cxx%*%t(A.th)
		ll0 = -sum(diag(SA.th%*%solve(S.th)))
        cnt = 0
		wA = which(A.th!=0)
		lA = length(wA)
		A = matrix(0,d*d,d*d)
		b = matrix(Cxy,d*d,1)
		lwi=0
		lwj=0
		for (id in 1:d){
   			wi = which(A.th[id,]!=0)
    		for (jd in 1:d){
        		cnt = cnt+1
        		A[cnt,lwi+(1:length(wi))] = Cxx[jd,wi]
        		A[cnt,lA+lwj+(1:(d-length(wi)))] = -S.th[-wi,jd]/2        
    		}
    		lwj = lwj+(d-length(wi))
    		lwi = lwi+length(wi)
		}
		A2.lasso = matrix(0,d,d)
		tmp = solve(A,b)[1:lA]
		lwi = 0
		for (id in 1:d){
    		wi = which(A.th[id,]!=0)
    		A2.lasso[id,wi] = tmp[lwi+(1:length(wi))]
    		lwi = lwi+length(wi)
		}        
		tmp = (m_1[,j]-(A2.lasso)%*%m[,j])/postmix[j] 	 
        tmp2 = Cyy + A2.lasso %*% Cxx %*% t(A2.lasso)-(A2.lasso %*%Cxy + t(A2.lasso%*%Cxy))
		S2.lasso = Cyy-Cxy%*%A2.lasso-A2.lasso%*%Cxy+A2.lasso%*%Cxx%*%t(A2.lasso)
		ll1 = -sum(diag(S2.lasso%*%solve(tmp2)))
		# print(paste("ll0 =",ll0, "ll1 =",ll1))
		
		# for (id in 1:d){
 			# w = which(theta$A[[j]][[1]][id,]!=0)
  			# Cxx_w = Cxx[w,w]
  			# Cxy_w = (Cxy)[w,id]
  			# A2.lasso[id,w] = t(Cxy_w)%*%solve(Cxx_w) 
  		# }
		tmp = (m_1[,j]-(A2.lasso)%*%m[,j])/postmix[j]
		tmp2 = Cyy + A2.lasso %*% Cxx %*% t(A2.lasso)-(A2.lasso %*%Cxy + t(A2.lasso%*%Cxy))
		A2[[j]][[1]] = A2.lasso
		if (sigma.diag) {
			sigma[[j]]=diag(diag(tmp2[1:d,1:d]),d)
		} else {
			w = which(abs(theta$sigma[[j]])>0)
			sigma[[j]] = matrix(0,d,d)
			sigma[[j]][w]=tmp2[1:d,1:d][w]
			#sigma[[j]]=tmp2[1:d,1:d]
		}
		moy[j,1:d] = tmp[1:d]
	}		
			
    if (p>0) {
	    list(A=A2,A0=moy,sigma=sigma,prior=prior,transmat=transmat)
	} 
}
