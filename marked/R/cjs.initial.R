#' Computes starting values for CJS p and Phi parameters
#' 
#' Computes starting values for Phi and p parameters from the
#' list of design matrices and the summarized data list including ch matrix and
#' first and last vectors. If any values are missing (NA) or abs(par)>5, they are
#' set to 0.
#' 
#' @param dml design matrix list for Phi and p
#' @param imat list containing chmat, first and last
#' @param link either "logit" (for cjs) or "probit" (for probitCJS)
#' @return list of initial parameter estimates
#' @author Jeff Laake
#' @keywords utility
cjs.initial=function(dml,imat,link="logit")
{
#   Create initial values for p using bernoulli glm with Manly-Parr approach
	message("Computing initial parameter estimates\n")
	num=nrow(imat$chmat)
	ind=matrix(c(1:num,imat$first+1,imat$last-1),ncol=3,nrow=num)
	ind=ind[ind[,2]<=ind[,3],]
	wts=unlist(apply(ind, 1, function(x, z) rep(z[x[1]], x[3]-x[2]+1),z = imat$freq))
	dep.values=unlist(apply(ind,1,function(x,z)z[x[1],x[2]:x[3]],z=imat$chmat))
	dd.indices=unlist(apply(ind,1,function(x) (x[1]-1)*(ncol(imat$chmat)-1)+x[2]:x[3]))-1
	x=dml[["p"]]$fe[dd.indices,]
    initial.p<-coef(glm.fit(x,dep.values,weights=wts,family=binomial(link=link)))
	initial.p[is.na(initial.p)]=0
	if(any(initial.p < -5 | initial.p > 5))initial.p=rep(0,length(initial.p))
	#   Create initial values for S using bernoulli glm assuming p=1
	last=imat$last+1
	last[last>ncol(imat$chmat)]=ncol(imat$chmat)
	ind=matrix(c(1:num,imat$first+1,last),ncol=3,nrow=num)
	ind=ind[ind[,2]<=ncol(imat$chmat),]
	wts=unlist(apply(ind, 1, function(x, z) rep(z[x[1]], x[3]-x[2]+1),z = imat$freq))
	dep.values=unlist(apply(ind,1,function(x,z)c(rep(1,(x[3]-x[2])),z[x[1],x[3]]),z=imat$chmat))
	dd.indices=unlist(apply(ind,1,function(x) (x[1]-1)*(ncol(imat$chmat)-1)+x[2]:x[3]))-1
	x=dml[["Phi"]]$fe[dd.indices,]
	initial.Phi<-coef(glm.fit(x,dep.values,weights=wts,family=binomial(link=link)))
	initial.Phi[is.na(initial.Phi)]=0
	if(any(initial.Phi < -5 | initial.Phi > 5))initial.Phi=rep(0,length(initial.Phi))
	return(list(Phi=initial.Phi,p=initial.p))
}
