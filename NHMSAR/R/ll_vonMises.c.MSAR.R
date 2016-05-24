ll_vonMises.c.MSAR <-
function(par,data,gamma,order,constr) {
	T = dim(data)[1]
	N.samples = dim(data)[2]
	d = dim(data)[3]
	f = 0
	res = deplie.c.VM(par,constr) 
	m = res$m
	kappa = res$kappa
	B = matrix(0,1,T-order)

	for (ex in 1:N.samples){
		para_comp = kappa[1]*exp(1i*m)
		if(order>0){
			for (o in 1:order) {
				para_comp = para_comp+kappa[o+1]*exp(1i*data[o:(T-order+o-1),ex,])
			}
		}	
		B = log_dens_Von_Mises(t(data[(order+1):T,ex,]),Arg(para_comp),abs(para_comp))
		# 
		f = f+sum(gamma[ex,]*B[1:length(gamma[ex,])])
		# 
	
	
}

f=-f;
return(f) 
}
