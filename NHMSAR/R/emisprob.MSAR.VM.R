emisprob.MSAR.VM <-
function(data,theta,covar=NULL) {
	
	d <- attributes(theta)$NbComp
	M <- attributes(theta)$NbRegimes
	order <- attributes(theta)$order
	label <- attributes(theta)$label
	data = as.matrix(data)
	T <- dim(as.matrix(data))[1]
	prob <- matrix(0,T-order,M)

	kappa = matrix(theta$kappa,M,order+1)
	mu=matrix(theta$mu,M,1)
#	if(d==1){
	for (j in 1:M) {
		para_comp = kappa[j,1]*exp(1i*mu[j])
		if(order>0){
			for (o in 1:order) {
				para_comp = para_comp+kappa[j,o+1]*exp(1i*data[o:(T-order+o-1),])
			}
		}
		prob[,j] = exp(log_dens_Von_Mises(t(data[(order+1):T,]),Arg(para_comp),abs(para_comp)))
		# 
   	}
    prob = t(prob)
    prob
}
