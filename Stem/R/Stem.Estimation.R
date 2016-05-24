`Stem.Estimation` <-
function(StemModel, precision=0.01, max.iter=50,flag.Gdiag=TRUE,flag.Sigmaetadiag=TRUE,cov.spat=Sigmastar.exp)	
{

z 		=  StemModel$data$z

p 		= StemModel$skeleton$p
n 		= StemModel$data$n
d 		= StemModel$data$d
r 		= StemModel$data$r

covariates 	= StemModel$data$covariates
covariates  = changedimension_covariates(covariates,d=d,r=r,n=n)
coordinates = StemModel$data$coordinates

phi_start 	= StemModel$skeleton$phi
phi_start$K = t(StemModel$skeleton$K)
n_par     	= length(unlist(phi_start))   

phi_start$logb 		= log(phi_start$sigma2eps/phi_start$sigma2omega)
phi_start$logtheta	= log(phi_start$theta)


####################
###EM algorithm while loop
####################
#output matrices (two columns are added for number of iteration and -2loglik)
parameters_mat   = matrix(0,nrow=max.iter,ncol=n_par+2)        
#colnames(parameters_mat) = c("n_iter","-2loglik",rep("k",dd*pp),"sigma2omega", "theta","logb",rep("beta",6),rep("G",pp*pp),rep("Sigmaeta",pp*pp),rep("m0",pp),rep("C0",pp*pp))
distance_mat    	= matrix(0,nrow=max.iter,ncol=1)
distancelog_mat	= c()
Q_mat       	= matrix(0,nrow=max.iter,ncol=2)   
iterNR 		= c()

converged_EM_1 	= FALSE
converged_EM_2 	= FALSE
n_iter_EM     	= 1

cat(paste("****************EM Algorithm - iteration n.",n_iter_EM),"\n")
while ((!converged_EM_1 | !converged_EM_2) && n_iter_EM < max.iter){
	step = kalman(	z            		= z,
			coordinates  	= coordinates,
         		p           		= p,
			n			= n,
			d			= d,
			r			= r,	
			phi_j        		= phi_start,
			max.iter     	= max.iter,
			precision       	= precision,
			covariates   	= covariates,
			Gdiag        		= flag.Gdiag,
			Sigmaetadiag 	= flag.Sigmaetadiag,
			cov.spat		= cov.spat
	)

	iterNR[n_iter_EM] 	= step$n_iter_NR
	step$phi$loglik 		= -2*step$phi$loglik 
	par           		= t(matrix(unlist(step$phi)))
	parameters_mat[n_iter_EM,] = cbind(n_iter_EM, par)

	if(n_iter_EM==1) {
  		prev_lik = 0
  		prev_par = rep(0,n_par)
   	} else {
  		prev_lik = (parameters_mat[n_iter_EM-1, 2])          #second column for -2loglik
  		prev_par = parameters_mat[n_iter_EM-1, -c(1,2)]   #no n_iter e -2loglik
  	}

  	if(n_iter_EM==1 | n_iter_EM==2) {
  		media_lik = 1
  	} else {
  		media_lik = mean(c(step$phi$loglik, prev_lik))
  	}

	dist_rel_num = sqrt(t(parameters_mat[n_iter_EM,-c(1,2)] - unlist(prev_par)) %*% (parameters_mat[n_iter_EM,-c(1,2)] - unlist(prev_par)))
	dist_rel_den = sqrt(t(unlist(prev_par)) %*% unlist(prev_par))
	dist_rel = dist_rel_num / dist_rel_den
	distance_mat[n_iter_EM,] = dist_rel

	diff_rel_loglik = abs(step$phi$loglik - prev_lik) / abs(prev_lik)
	distancelog_mat[n_iter_EM] = diff_rel_loglik

	###Check the convergence!
	converged_EM_1 = diff_rel_loglik < precision
	converged_EM_2 = dist_rel < precision

	Q_mat[n_iter_EM,] = c(unlist(step$Q_prev),unlist(step$Q_new))

	###Updating the iteration number and the parameter vector
	phi_start 	= step$phi[-1]
	n_iter_EM 	= n_iter_EM + 1
	cat(paste("****************EM Algorithm - iteration n.",n_iter_EM),"\n")

} # here the while loop ends

phi_start$theta = exp(phi_start$logtheta)  
phi_start$sigma2eps = exp(phi_start$logb) * phi_start$sigma2omega  
phi_start = phi_start[-which(names(phi_start) == "logtheta")]
phi_start = phi_start[-which(names(phi_start) == "K")]
phi_start = phi_start[-which(names(phi_start) == "logb")]
phi.estimated= phi_start

StemModel$estimates$phi.hat = phi.estimated 
StemModel$estimates$y.smoothed = step$m.smoother 
StemModel$estimates$loglik = (parameters_mat[(n_iter_EM-1),2])*(-2) 
convergence.par 			= list(conv.log = converged_EM_1,
						conv.par = converged_EM_2,
						iterEM   = n_iter_EM-1,
						iterNR   = iterNR)
StemModel$estimates$convergence.par = convergence.par
return (StemModel)		
}

