#' Likelihood function for Jolly-Seber model using Schwarz-Arnason POPAN
#' formulation
#' 
#' For a given set of parameters and data, it computes -2*log Likelihood value
#' but does not include data factorials. Factorials for unmarked are not needed
#' but are included in final result by \code{\link{js}} so the result matches
#' output from MARK for the POPAN model.
#' 
#' This functions uses \code{\link{cjs.lnl}} and then supplements with the
#' remaining calculations to compute the likelihood for the POPAN formulation
#' (Arnason and Schwarz 1996) of the Jolly-Seber model.
#' 
#' @param par vector of parameter values
#' @param model_data a list that contains: 1)imat-list of vectors and matrices constructed by
#' \code{\link{process.ch}} from the capture history data, 2)Phi.dm design matrix for Phi constructed by \code{\link{create.dm}},
#' 3)p.dm design matrix for p constructed by \code{\link{create.dm}}, 4)pent.dm design matrix for probability of entry constructed by \code{\link{create.dm}},
#' 5) N.dm design matrix for estimates of number of animals not caught from
#' super-population constructed by \code{\link{create.dm}},
#' 6)Phi.fixed matrix with 3 columns: ch number(i), occasion number(j),
#' fixed value(f) to fix phi(i,j)=f, 7) p.fixed matrix with 3 columns: ch number(i), occasion number(j), 
#' 8) pent.fixed matrix with 3 columns: ch number(i), occasion number(j), fixed value(f) to fix pent(i,j)=f, and
#' 9) time.intervals intervals of time between occasions if not all 1
#' fixed value(f) to fix p(i,j)=f
#' @param debug if TRUE will printout values of \code{par} and function value
#' @param nobstot number of unique caught at least once by group if applicable
#' @param jsenv environment for js to update iteration counter
#' @return -log likelihood value, excluding data (ui) factorials which are added in js after optimization to match MARK
#' @author Jeff Laake
#' @references Schwarz, C. J., and A. N. Arnason. 1996. A general methodology
#' for the analysis of capture-recapture experiments in open populations.
#' Biometrics 52:860-873.
js.lnl=function(par,model_data,debug=FALSE,nobstot,jsenv)
{
	get.pent=function(beta,dm,nocc)
	{
		pents=cbind(rep(1,nrow(dm)/(nocc-1)),exp(matrix(as.vector(dm%*%beta),ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE)))
		pents=pents/apply(pents,1,sum)
		return(pents)
	}
	get.p=function(beta,dm,nocc)
	{
		ps=cbind(rep(1,nrow(dm)/(nocc-1)),plogis(matrix(as.vector(dm%*%beta),ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE)))
		return(ps)
	}
# compute Phi matrix from parameters (beta) and list of design matrices (dm)
# created by function create.dm
	get.Phi=function(beta,dm,nocc)
	{
		Phis=cbind(rep(1,nrow(dm)/(nocc-1)),plogis(matrix(as.vector(dm%*%beta),ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE)))
		return(Phis)
	}
	
if(debug)cat("par = ",par,"\n")
f_eval=get("markedfunc_eval",envir=jsenv)+1
assign("markedfunc_eval", f_eval, envir = jsenv)
# initialize constants and parameter vectors
nocc=model_data$imat$nocc
nphi=ncol(model_data$Phi.dm)
np=ncol(model_data$p.dm)
npent=ncol(model_data$pent.dm)
nN=ncol(model_data$N.dm)
beta.phi=par[1:nphi]
beta.p=par[(nphi+1):(nphi+np)]
beta.pent=par[(nphi+np+1):(nphi+np+npent)]
beta.N=par[(nphi+np+npent+1):(nphi+np+npent+nN)]
# create Phi and p beta matrices excluding first occasion on p
Phibeta=matrix(as.vector(model_data$Phi.dm%*%beta.phi),ncol=nocc-1,nrow=nrow(model_data$Phi.dm)/(nocc-1),byrow=TRUE)
pbeta=matrix(as.vector(model_data$p.dm%*%beta.p),ncol=nocc,nrow=nrow(model_data$p.dm)/(nocc),byrow=TRUE)
if(is.null(model_data$time.intervals))model_data$time.intervals=rep(1,nocc-1)
# compute CJS portion of the likelihood
cjslnl=.Fortran("cjs",as.double(model_data$imat$chmat),as.double(Phibeta),as.double(pbeta[,-1]),
           as.double(model_data$imat$first),as.double(model_data$imat$last),as.double(model_data$imat$freq),
           as.integer(model_data$imat$loc),as.double(model_data$Phi.fixed),as.double(model_data$p.fixed[model_data$p.fixed[,2]!=1,,drop=FALSE]),
           as.double(model_data$time.intervals),as.integer(nrow(model_data$imat$chmat)),           
           as.integer(ncol(model_data$imat$chmat)),as.integer(nrow(model_data$Phi.fixed)),
           as.integer(nrow(model_data$p.fixed[model_data$p.fixed[,2]!=1,,drop=FALSE])),lnl=double(1),
           p0=double(nrow(model_data$imat$chmat)),PACKAGE="marked")
# next add on likelihood component for first capture
pents=get.pent(beta.pent,model_data$pent.dm,nocc)
pents.dummy=pents[model_data$imat$freq==0,]
ps=plogis(matrix(as.vector(model_data$p.dm%*%beta.p),ncol=nocc,nrow=nrow(model_data$p.dm)/(nocc),byrow=TRUE))
ps.dummy=ps[model_data$imat$freq==0,]
Phis=get.Phi(beta.phi,model_data$Phi.dm,nocc)
p.occ=ps[cbind(1:nrow(pents),model_data$imat$first)]
ps[,nocc]=0
Phis=cbind(Phis[,2:nocc],rep(1,nrow(Phis)))
entry.p=(1-ps)*Phis*(1-model_data$imat$First)+model_data$imat$First
# Looping faster because k<<<n
#Estar=t(apply(entry.p,1,function(x) rev(cumprod(rev(x)))))
Estar=matrix(0,ncol=nocc,nrow=nrow(ps))
Estar[,nocc]=entry.p[,nocc]
for(j in (nocc-1):1)
	Estar[,j]=Estar[,j+1]*entry.p[,j]
entry.p=rowSums(Estar*pents*(1-model_data$imat$Fplus))*p.occ
lnl=cjslnl$lnl-sum(model_data$imat$freq*log(entry.p))
# next add on likelihood component for those not caught from dummy 1000000,0100000,...data 
# return complete likelihood value except that calling function js adds the ui factorials to match
# POPAN output from MARK
Ns=exp(as.vector(model_data$N.dm%*%beta.N))
for (i in 1:length(Ns))
{
  index0=(i-1)*nocc+1
  index1=i*nocc
  ps=ps.dummy[index0:index1,]
  pents=pents.dummy[index0:index1,]
  lnl=lnl-Ns[i]*log(sum(diag(pents)*cjslnl$p0[model_data$imat$freq==0][index0:index1]*(1-diag(ps))))
  lnl=lnl-lfactorial(nobstot[i]+Ns[i])+lfactorial(Ns[i])
}
if(debug)
{
	cat("-2lnl = ",2*lnl,"\n")
} else
    if((f_eval-100*floor(f_eval/100))==0)
    {
		cat("\r Number of evaluations: ",f_eval," -2lnl:",formatC(2*lnl,digits=10))
		flush.console()
    }	
return(lnl)
}
