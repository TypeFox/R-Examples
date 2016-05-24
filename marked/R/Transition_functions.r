#' HMM Transition matrix functions
#' 
#' Functions that compute the transition matrix for various models. Currently only CJS and MS models
#' are included.
#'  
#' @param pars list of real parameter values for each type of parameter
#' @param m number of states
#' @param F initial occasion vector 
#' @param T number of occasions
#' @aliases cjs_gamma ms_gamma ms2_gamma
#' @export cjs_gamma ms_gamma ms2_gamma
#' @return array of id and occasion-specific transition matrices - Gamma in Zucchini and MacDonald (2009)
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. 275p. 
cjs_gamma=function(pars,m,F,T) 
 {
	phimat=array(0,c(nrow(pars$Phi),T-1,m,m))
 	# create 4-d array with a matrix for each id and occasion
 	# from pars$Phi which is a matrix of id by occasion survival probabilities 	
 	value=.Fortran("cjsgam",as.double(pars$Phi),as.integer(nrow(pars$Phi)),
			as.integer(F),as.integer(T),phimat=double(nrow(pars$Phi)*(T-1)*4),PACKAGE="marked")
	dim(value$phimat)=c(nrow(pars$Phi),T-1,2,2)
	value$phimat
}
cjs1tl_gamma=function(pars,m,F,T) 
{
	# create 4-d array with a matrix for each id and occasion
	# from pars$Phi which is a matrix of id by occasion survival probabilities 	
	tmat=array(0,c(nrow(pars$Phi),T-1,m,m))
	value=.Fortran("cjs1tlgam",as.double(pars$Phi),as.double(pars$tau),as.integer(nrow(pars$Phi)),
			as.integer(F),as.integer(T),tmat=double(nrow(pars$Phi)*(T-1)*m^2),PACKAGE="marked")
	dim(value$tmat)=c(nrow(pars$Phi),T-1,m,m)
	value$tmat
}
cjs2tl_gamma=function(pars,m,F,T) 
{
	  # create 4-d array with a matrix for each id and occasion
	  # from pars$Phi which is a matrix of id by occasion survival probabilities 	
	  tmat=array(0,c(nrow(pars$Phi),T-1,m,m))
	  value=.Fortran("cjs2tlgam",as.double(pars$Phi),as.double(pars$tau),as.integer(nrow(pars$Phi)),
			  as.integer(F),as.integer(T),tmat=double(nrow(pars$Phi)*(T-1)*m^2),PACKAGE="marked")
	  dim(value$tmat)=c(nrow(pars$Phi),T-1,m,m)
	  value$tmat
}
ms_gamma=function(pars,m,F,T) 
{
	 # create 4-d array with a matrix for each id and occasion
	 # from pars$S  which is a matrix of id by occasion survival probabilities 
	 # and pars$Psi which is state transitions
	if(is.list(m))m=m$ns*m$na+1
	value=.Fortran("msgam",as.double(pars$S),as.double(pars$Psi),as.integer(nrow(pars$S)),as.integer(m),
			 as.integer(F),as.integer(T),tmat=double(nrow(pars$S)*(T-1)*m^2),PACKAGE="marked")
	 dim(value$tmat)=c(nrow(pars$S),T-1,m,m)
	 value$tmat
 }
 ms2_gamma=function(pars,m,F,T) 
 {
	 ns=m$ns*m$na+1
	 # create 4-d array with a matrix for each id and occasion
	 # from pars$Phi which is a matrix of id by occasion survival probabilities 	
	 value=.Fortran("ms2gam",as.double(pars$S),as.double(pars$Psi),as.double(pars$alpha),as.integer(nrow(pars$S)),
			 as.integer(ns),as.integer(m$na),as.integer(m$ns),as.integer(F),as.integer(T),
			 tmat=double(nrow(pars$S)*(T-1)*ns^2),PACKAGE="marked")
	 dim(value$tmat)=c(nrow(pars$S),T-1,ns,ns)
	 value$tmat
 }
 mvms_gamma=function(pars,m,F,T) 
 {
	 # create 4-d array with a matrix for each id and occasion
#	 # from pars$Phi which is a matrix of id by occasion survival probabilities 
	 # and pars$Psi which is state transitions
	 value=.Fortran("msgam",as.double(pars$Phi),as.double(pars$Psi),as.integer(nrow(pars$Phi)),as.integer(m),
			 as.integer(F),as.integer(T),tmat=double(nrow(pars$Phi)*(T-1)*m^2),PACKAGE="marked")
	 dim(value$tmat)=c(nrow(pars$Phi),T-1,m,m)
	 value$tmat
 }
# R versions of the FORTRAN code
#  cjs_gamma=function(pars,m,F,T) 
#  {
#  	# create 4-d array with a matrix for each id and occasion
#  	# from pars$Phi which is a matrix of id by occasion survival probabilities 	
#  	phimat=array(NA,c(nrow(pars$Phi),T-1,m,m))
#  	for (i in 1:nrow(phimat))
#  		for(j in F[i]:(T-1))
#  		{
#  			phi=pars$Phi[i,j]
#  			phimat[i,j,,]=matrix(c(phi,1-phi,0,1),nrow=2,ncol=2,byrow=TRUE)
#  		}
#  	phimat
#  }
#ms_gamma=function(pars,m,F,T) 
#{
#	# create an 4-d array with a matrix for each id and occasion for S from pars$S 
#	# which is a matrix of id by occasion x state survival probabilities
#	if(is.list(m))m=m$ns*m$na+1
#	phimat=array(NA,c(nrow(pars$S),T-1,m,m))
#	for (i in 1:nrow(phimat))
#	{
#		for(j in F[i]:(T-1))
#		{
#			s=pars$S[i,((j-1)*(m-1)+1):(j*(m-1))]
#			smat=matrix(s,ncol=length(s),nrow=length(s))
#			phimat[i,j,,]=rbind(cbind(smat,1-s),c(rep(0,length(s)),1))
#		}
#	}
#	# create a 4-d array from pars$Psi which is a matrix of id by occasion x state^2 
#	# non-normalized Psi probabilities which are normalized to sum to 1. 			
#	psimat=array(NA,c(nrow(pars$Psi),T-1,m,m))
#	for (i in 1:nrow(psimat))
#	{
#		for(j in F[i]:(T-1))
#		{
#			psi=pars$Psi[i,((j-1)*(m-1)^2+1):(j*(m-1)^2)]
#           psix=matrix(psi,ncol=sqrt(length(psi)),byrow=TRUE)
#			psix=psix/rowSums(psix)
#			psimat[i,j,,]=rbind(cbind(psix,rep(1,nrow(psix))),rep(1,nrow(psix)+1))
#		}
#	}
#	# The 4-d arrays are multiplied and returned
#	phimat*psimat
#}
# Following wasn't really needed since it only made change from S to Phi
# and made sure it only used F[i]<T.The latter was put into FORTRAN code for msgamma 
#mvms_gamma=function(pars,m,F,T) 
#{
#	# create an 4-d array with a matrix for each id and occasion for S from pars$S 
#	# which is a matrix of id by occasion x state survival probabilities
#	if(is.list(m))m=m$ns*m$na+1
#	phimat=array(0,c(nrow(pars$Phi),T-1,m,m))
#	for (i in 1:nrow(phimat))
#	{
#		if(F[i]<=(T-1))
#		for(j in F[i]:(T-1))
#		{
#			s=pars$Phi[i,((j-1)*(m-1)+1):(j*(m-1))]
#			smat=matrix(s,ncol=length(s),nrow=length(s))
#			phimat[i,j,,]=rbind(cbind(smat,1-s),c(rep(0,length(s)),1))
#		}
#	}
#	# create a 4-d array from pars$Psi which is a matrix of id by occasion x state^2 
#	# non-normalized Psi probabilities which are normalized to sum to 1. 
#	psimat=array(0,c(nrow(pars$Psi),T-1,m,m))
#	for (i in 1:nrow(psimat))
#	{
#		if(F[i]<=(T-1))
#		for(j in F[i]:(T-1))
#		{
#			psi=pars$Psi[i,((j-1)*(m-1)^2+1):(j*(m-1)^2)]
#           psix=matrix(psi,ncol=sqrt(length(psi)),byrow=TRUE)
#			psix=psix/rowSums(psix)
#			psimat[i,j,,]=rbind(cbind(psix,rep(1,nrow(psix))),rep(1,nrow(psix)+1))
#		}
#	}
#	# The 4-d arrays are multiplied and returned
#	phimat*psimat
#}
#ms2_gamma=function(pars,m,F,T) 
#{
#	# create an 4-d array with a matrix for each id and occasion for S from pars$S 
#	# which is a matrix of id by occasion x state survival probabilities
#	ns=m$ns*m$na+1
#	phimat=array(NA,c(nrow(pars$S),T-1,ns,ns))
#	for (i in 1:nrow(phimat))
#	{
#		for(j in F[i]:(T-1))
#		{
#			s=pars$S[i,((j-1)*(ns-1)+1):(j*(ns-1))]
#			smat=matrix(s,ncol=length(s),nrow=length(s))
#			phimat[i,j,,]=rbind(cbind(smat,1-s),c(rep(0,length(s)),1))
#		}
#	}
#	# create a 4-d array from pars$Psi which is a matrix of id by occasion x state^2 
#	# non-normalized Psi probabilities which are normalized to sum to 1. 			
#	psimat=array(NA,c(nrow(pars$Psi),T-1,ns,ns))
#	for (i in 1:nrow(psimat))
#	{
#		for(j in F[i]:(T-1))
#		{
#			psi=pars$Psi[i,((j-1)*m$ns^2+1):(j*m$ns^2)]
#			psix=matrix(psi,ncol=sqrt(length(psi)),byrow=TRUE)
#			psix=psix/rowSums(psix)
#			psix=matrix(rep(apply(psix,2,rep,m$na),m$na),nrow=m$ns*m$na)
#			alpha=pars$alpha[i,((j-1)*m$na^2+1):(j*m$na^2)]
#			alphax=matrix(alpha,ncol=sqrt(length(alpha)),byrow=TRUE)
#			alphax=alphax/rowSums(alphax)
#			alphax=matrix(apply(alphax,2,function(x) rep(rep(x,each=m$ns),times=m$ns)),nrow=m$ns*m$na) 
#			psimat[i,j,,]=rbind(cbind(psix*alphax,rep(1,nrow(psix))),rep(1,nrow(psix)+1))
#		}
#	}
#	# The 4-d arrays are multiplied and returned
#	phimat*psimat
#}
#cjs2tl_gamma=function(pars,m,F,T) 
#{
#	# create 4-d array with a matrix for each id and occasion
#	# from pars$Phi which is a matrix of id by occasion survival probabilities 	
#	tmat=array(0,c(nrow(pars$Phi),T-1,m,m))
##     normalize to sum to 1
#	for (i in 1:nrow(tmat))
#	{
#		tau=matrix(pars$tau[i,],ncol=4,byrow=TRUE)
#		tau=tau/rowSums(tau)
#		tau0=matrix(0,ncol=2,nrow=nrow(tau))   
#		tau0[,1]=ifelse((tau[,2]+tau[,4])>0,tau[,4]/(tau[,2]+tau[,4]),0)  
#		tau0[,2]=ifelse((tau[,3]+tau[,4])>0,tau[,4]/(tau[,3]+tau[,4]),0)  
#		for(j in F[i]:(T-1))
#		{
#			phi=pars$Phi[i,j]
#			tmat[i,j,,]=matrix(c(phi*tau[j,],1-phi,
#							0,phi*(1-tau0[j,1]),0,phi*tau0[j,1],1-phi,
#							0,0,phi*(1-tau0[j,2]),phi*tau0[j,2],1-phi,								 
#							0,0,0,phi,1-phi,
#							0,0,0,0,1),nrow=5,ncol=5,byrow=TRUE)
#		}
#	}
#	tmat
#}
#cjs1tl_gamma=function(pars,m,F,T) 
#{
#	# create 4-d array with a matrix for each id and occasion
#	# from pars$Phi which is a matrix of id by occasion survival probabilities 	
#	tmat=array(0,c(nrow(pars$Phi),T-1,m,m))
#	for (i in 1:nrow(tmat))
#		for(j in F[i]:(T-1))
#		{
#			phi=pars$Phi[i,j]
#			tau=pars$tau[i,j]
#			tmat[i,j,,]=matrix(c(phi*(1-tau),tau*phi,1-phi,0,phi,1-phi,0,0,1),nrow=3,ncol=3,byrow=TRUE)
#		}
#	tmat
#}
