Mstep.classif <-
function(data,C,order)  {  
	T=dim(data)[1]
	N.samples = dim(as.array(data))[2] 
	d = dim(as.array(data))[3]
	if(is.null(d)|is.na(d)) {d = 1}
    M <- length(unique(c(C)))
    p <- order
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
        
    #T = length(o:(dim(data)[1]-order+o))
    tmp = order:dim(C)[1]
    C = matrix(C[tmp,,],length(tmp),N.samples)
    T = dim(C)[1]
#	exp_num_trans = 0
#	exp_num_visit = 0
#	exp_num_visits1 = 0
	postmix = 0
	m = matrix(0,d*order,M) ; m_1 = m ; c = matrix(0,M,1) ; s=c ; 
    op = array(0,c(d*order,d*order,M) ); op_1 = op ; op_2 = op_1 ;

	for (ex in 1:N.samples) {
		obs = array(data2[,,ex],c(d*order,T))
		gamma = matrix(0.001,T,M)
		for (km in 1:M) {gamma[which(C[,ex]==km),km] = 0.999}

		postmix = postmix+apply(gamma,2,sum)
		obs = t(obs) 
		    
		for (j in 1:M) {
			w = matrix(gamma[1:(T-1),j], 1,T-1)
		    wobs = obs[1:(T-1),] * repmat(t(w),1,d*order)
            wobs_1 = obs[2:T,] * repmat(t(w),1,d*order)
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
     
     moy <- array(0,c(M,d))
     Sigma <- array(0,c(d,d,M))
	 A2 <-array(0,c(M,order,d^2))
																		
	for (j in 1:M) { 
		Cxx = postmix[j]*op[,,j] - m[,j]%*%t(m[,j])
		if (det(Cxx)<=0) {Cxx = Cxx+1e-4*diag(1,d)}
        Cxy = postmix[j]*op_1[,,j] - m[,j]%*%t(m_1[,j])
        #A2tmp = t(Cxy)%*%solve(Cxx) 
        A2tmp = t(solve(Cxx,Cxy)) 
        tmp = (m_1[,j]-(A2tmp)%*%m[,j])/postmix[j] 
        moy[j,1:d] = tmp[1:d] 
        op_1[,,j] = t(as.matrix(op_1[,,j]))
        tmp2 = (op_2[, , j] + (A2tmp) %*% op[, , j] %*% t(A2tmp) - 
                   ((A2tmp) %*% t(op_1[, , j]) + t((A2tmp) %*% t(op_1[, , j]))))/postmix[j] - tmp %*% t(tmp)
        Sigma[1:d,1:d,j]=tmp2[1:d,1:d]
        for(kp in 1:order){
            A2[j,kp,] = A2tmp[1:d,((kp-1)*d+1):(kp*d)]
        }
	}
		
    if(d>1){
        Atmp=list()
        sigma=list()
        for(j in 1:M){
        	Atmp[[j]]=list()
        	for(kp in 1:order){
        		Atmp[[j]][[kp]]=matrix(A2[j,kp,],d,d)
        	}
        	sigma[[j]]=Sigma[,,j]
        }
     }else{
		Atmp=matrix(A2,M,order)
		sigma=matrix(Sigma,byrow=TRUE,nrow=M)
	}
	LLy = -T*N.samples*d/2-T*d/2 # observations
	LLx = 0 # classification
	for (i in 1:M){
		ni = sum(C==i)
		LLy = LLy-ni*log(det(as.matrix(sigma[[i]])))
		for (j in 1:M) {
			nij = sum(C[1:(T-1),]==i & C[2:(T),]==j)
			if (nij>0) {LLx = LLx + sum(nij)*log(nij/ni)}
		}
	}
	if (M>1) {LL = LLy+LLx}
	else {LL = LLy}
	list(A=Atmp,A0=moy,sigma=sigma,LL=LL)
	
}
