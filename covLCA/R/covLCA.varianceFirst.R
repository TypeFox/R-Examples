covLCA.varianceFirst <-
function(prior,probs,rgivy,R,S1,J,K.j,S2,x,y,z)
{
	library(Matrix)
	dd.bet <- covLCA.dQdBeta(rgivy,prior,x) #A: a list containing the gradient and the hessian of Q_beta wrt betas, beta_jp
	
	hessAlist=vector(length=J,mode="list")
	hessGlist=vector(length=J,mode="list")
	hessGAlist=vector(length=J,mode="list")
	hessAGlist=vector(length=J,mode="list")
	
	for (m in 1:J) #A: for each manifest variable (i.e., for each subset of parameters alphas and gammas)
	{
		dd.gam <- covLCA.dQdGamma(rgivy,probs,y,K.j,m) #A: a list containing the gradient and the hessian of Q_omega wrt gammas, gamma_(m)jk
		dd.alph <- covLCA.dQdAlpha(rgivy,probs,z,K.j,m,y,S2) #A: a list containing the gradient and the hessian of Q_omega wrt alphas, alpha_(m)qk
		dd.alph.gam <- covLCA.dQdAlphaGamma(rgivy,probs,z,K.j,m,S2)
							
		hessGlist[[m]]=dd.gam$hess #A: a list where element m is the hessian matrix of Q_omega_m wrt gamma_m
		hessAlist[[m]]=dd.alph$hess #A: a list where element m is the hessian matrix of Q_omega_m wrt alpha_m
		hessAGlist[[m]]=t(dd.alph.gam) #A: a list where element m is the hessian matrix of Q_omega_m wrt alpha_m and gamma_m (transpose of element of hessGAlist)
		hessGAlist[[m]]=dd.alph.gam #A: a list where element m is the hessian matrix of Q_omega_m wrt gamma_m and alpha_m
	}
	
	
	GG=as.matrix(bdiag(hessGlist)) ; AA=as.matrix(bdiag(hessAlist)) ; AG=as.matrix(bdiag(hessAGlist)) ; GA=as.matrix(bdiag(hessGAlist))
	zero1=matrix(0,nrow=(R-1)*S1,ncol=J*R*(K.j[1]-1)+J*S2*(K.j[1]-1)) #Does not need to be adapted
	
	EIc=rbind(cbind(dd.bet$hess,zero1),cbind(matrix(0,nrow=J*R*(K.j[1]-1),ncol=(R-1)*S1),GG,GA),cbind(matrix(0,nrow=J*S2*(K.j[1]-1),ncol=(R-1)*S1),AG,AA)) #Does not need to be adapted to "S2=1"
	
	return(EIc)
}
