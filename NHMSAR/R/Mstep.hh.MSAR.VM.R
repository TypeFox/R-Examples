Mstep.hh.MSAR.VM <-
function(data,theta,FB,constr=0)  {  
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
	}
      
     # Markov chain
     prior = normalise(exp_num_visits1)
     prior=matrix(prior,M,1)
     transmat = mk_stochastic(exp_num_trans)
     moy = matrix(0,M,1)
     kappa = matrix(0,M,p+1)
 
     if (!is.complex(theta$kappa)){     
     	for (j in 1:M) {
			a  = plie.VM(theta$mu[j],theta$kappa[j,])
			resopt = ucminf(a,fn=ll_vonMises.MSAR,data=data,gamma=FB$probS[,,j],order=p,gr = NULL, control = list(trace=FALSE)) 
			a.dep = deplie.VM(resopt$par)
			moy[j] =  a.dep$m
			kappa[j,] = a.dep$kappa
		
      	}
      } else {
      	for (j in 1:M) {
      		a  = plie.c.VM(theta$mu[j],theta$kappa[j,],constr)
			resopt = ucminf(a,fn=ll_vonMises.c.MSAR,data=data,gamma=FB$probS[,,j],order=p,constr=constr,gr = NULL, control = list(trace=FALSE))
			 
			a.dep = deplie.c.VM(resopt$par,constr)
			moy[j] =  a.dep$m
			kappa[j,] = a.dep$kappa  
				  
		}
      }
 
        if (min(postmix)<1e-6) {stop("error : smoothing probabilities are to small, in one regime at least. You should revise initialisation.")}															
	   list(mu=moy,kappa=kappa,prior=prior,transmat=transmat)
}
