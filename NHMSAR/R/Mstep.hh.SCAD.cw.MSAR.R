
# component wise optimisation (ncvreg)
Mstep.hh.SCAD.cw.MSAR <-
function(data,theta,FB,lambda1=.1,lambda2=.1,penalty="SCAD",par=NULL)  {  


	T=dim(data)[1]
	N.samples = dim(as.array(data))[2] 
	d = dim(as.array(data))[3]
	if(is.null(d)|is.na(d)) {d = 1}
    M <- attributes(theta)$NbRegimes
    p <- attributes(theta)$order
	order <- max(p,1)
	if (length(lambda1)==1) {lambda1 = matrix(lambda1,1,M)}
	if (length(lambda2)==1) {lambda2 = matrix(lambda2,1,M)}
	data = array(data,c(T,N.samples,d))
	data2 = array(0,c(order*d,T-order+1,N.samples))
	cpt=1
    for (o in order:1) {
    	for (kd in 1:d) {
    		 data2[cpt,,] = data[o:(T-order+o),,kd]
    		 cpt =cpt+1
    	}
    }
        #browser()
    T = length(o:(dim(data)[1]-order+o))
    N = T*N.samples
	exp_num_trans = 0
	exp_num_visit = 0
	exp_num_visits1 = 0
	postmix = 0
	m = matrix(0,d*order,M) ; m_1 = m ; c = matrix(0,M,1) ; s=c ; 
    op = array(0,c(d*order,d*order,M) ); op_1 = op ; op_2 = op_1 ;
    wx = array(0,c(N.samples*(T-1),d,M)) 
    wy = array(0,c(N.samples*(T-1),d,M))

    Cxx = array(0,c(d*order,d*order,M) )
    Cxy = array(0,c(d*order,d*order,M) )
    Cyy = array(0,c(d*order,d*order,M) )
	for (ex in 1:N.samples) {
		obs = t(array(data2[,,ex],c(d*order,T)))
		xit  = array(FB$probSS[,,,ex],c(M,M,T-2))
		gamma = matrix(FB$probS[ex,1:(T-1),],T-1,M)
		exp_num_trans = exp_num_trans+apply(xit,c(1,2),sum)
		exp_num_visits1 = exp_num_visits1+gamma[1,]
		postmix = postmix+apply(gamma,2,sum)
		for (j in 1:M) { # attention, ici on n'a qu'un ex
		   	w = matrix(gamma[,j], 1,T-1)
		    wx[(ex-1)*(T-1)+(1:(T-1)),,j] = obs[1:(T-1),] * repmat(t(w),1,d*order)
            wy[(ex-1)*(T-1)+(1:(T-1)),,j]= obs[2:T,] * repmat(t(w),1,d*order)
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
	 Areg=list()
	 Sreg = list()	
	 sigma.inv=list()
	 a = 3.7
	 for (j in 1:M) { # attention, ici on n'a qu'un ex
		Cxx[,,j] = postmix[j]*op[,,j] - m[,j]%*%t(m[,j])
      	Cxy[,,j] = postmix[j]*op_1[,,j] - m[,j]%*%t(m_1[,j])
      	Cyy[,,j] = postmix[j]*op_2[,,j] - m_1[,j]%*%t(m_1[,j])
      	A2tmp = t(Cxy[,,j])%*%solve(Cxx[,,j]) 
		A = theta$A[[j]][[1]]
		S = theta$sigma[[j]]		  
		S.emp = (Cyy[,,j]-Cxy[,,j]%*%t(A)-A%*%t(Cxy[,,j])+A%*%Cxx[,,j]%*%t(A))/postmix[j]^2
		res.gl = list()
		#S = S.emp
		#res.gl$wi = solve(S)
		if (!is.null(par$sigma.inv)) {res.gl$wi = par$sigma.inv[[j]]
		} else {res.gl$wi = solve(S)}
      	w = matrix(0,d,d)
		if (lambda1[j]>0 & penalty=="SCAD"){ 
		 	  w.tmp = apply(matrix(c(a*lambda1[j]-abs(res.gl$wi),rep(0,d*d)),d*d,2),1,max)
		 	  w = (lambda1[j]*(as.numeric(abs(res.gl$wi)<=lambda1[j])+w.tmp/((a-1)*lambda1[j])*as.numeric(abs(res.gl$wi)>lambda1[j])))
	    } else {w = matrix(lambda1[j],d,d)}
		
		if (lambda2[j]>0 & penalty=="SCAD") {
			  A.scad = matrix(0,d,d)
			  for (id in 1:d){
  				  A.scad[id,] = ncvreg(wx[,,j],wy[,id,j],family = "gaussian",penalty="SCAD",
  				                     lambda=lambda2[j])$beta[2:(d+1)]
  		}
  		A = A.scad 
		} else if (lambda2[j]>0 &  penalty=="ridge") {
        	omega=matrix(lambda2[j],d,d) 
			omega = omega-diag(diag(omega))
			res=ucminf(c(A), fn=ll.pen2,gr=grad.ll.pen2,S=S,Cxx=Cxx[,,j],Cxy=Cxy[,,j],Cyy=Cyy[,,j],omega=omega)
			A = res$par
		} else {A  = A2tmp}        
		Areg[[j]] = list()
		Areg[[j]][[1]] = A
		if (lambda1[j]>0) {	
			w = N*w
			res.gl = glasso(S, rho=matrix(w,d,d),penalize.diagonal=FALSE,maxit=10)
			sigma.inv[[j]] = res.gl$wi 
			Sreg[[j]] = res.gl$w
	  }
		else {
			tmp = (m_1[,j]-(A)%*%m[,j])/postmix[j] 
      op_1[,,j] = t(as.matrix(op_1[,,j]))
      Sreg[[j]] = (op_2[, , j] + (A) %*% op[, , j] %*% t(A) - 
                   ((A) %*% t(op_1[, , j]) + t((A) %*% t(op_1[, , j]))))/postmix[j] - tmp %*% t(tmp)
                   sigma.inv[[j]] = solve(Sreg[[j]])
		}
		moy[j,] = (m_1[,j]-A%*%m[,j])/postmix[j]
		
	}		
																
  if (p>0) {
	    list(A=Areg,A0=moy,sigma=Sreg,prior=prior,transmat=transmat,sigma.inv=sigma.inv)
	} else {
		list(sigma=Sreg,A0=moy,prior=prior,transmat=transmat,sigma.inv=sigma.inv)
	}
}
