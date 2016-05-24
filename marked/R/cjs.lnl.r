#' Likelihood function for Cormack-Jolly-Seber model
#' 
#' For a given set of parameters and data, it computes -log Likelihood value.
#' 
#' This function uses a FORTRAN subroutine (cjs.f) to speed up computation of the likelihood but the
#' result can also be obtained wholly in R with a small loss in precision. See R code below.
#' The R and FORTRAN code uses the likelihood formulation of Pledger et al.(2003).
#' \preformatted{
#' get.p=function(beta,dm,nocc,Fplus) 
#' {
#' # compute p matrix from parameters (beta) and list of design matrices (dm) 
#' # created by function create.dm
#'   ps=cbind(rep(1,nrow(dm)/(nocc-1)),
#'      matrix(dm%*%beta,ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE))
#'   ps[Fplus==1]=plogis(ps[Fplus==1]) 
#'   return(ps) 
#' }
#' get.Phi=function(beta,dm,nocc,Fplus) 
#' {
#' # compute Phi matrix from parameters (beta) and list of design matrices (dm)
#' # created by function create.dm 
#'   Phis=cbind(rep(1,nrow(dm)/(nocc-1)),
#'      matrix(dm%*%beta,ncol=nocc-1,nrow=nrow(dm)/(nocc-1),byrow=TRUE))
#'   Phis[Fplus==1]=plogis(Phis[Fplus==1]) 
#'   return(Phis) 
#' }
#'#################################################################################
#'# cjs.lnl - computes likelihood for CJS using Pledger et al (2003)
#'# formulation for the likelihood. This code does not cope with fixed parameters or 
#'# loss on capture but could be modified to do so. Also, to work directly with cjs.r and
#'# cjs.accumulate call to process.ch would have to have all=TRUE to get Fplus and L.
#'# Arguments: 
#'# par             - vector of beta parameters 
#'# imat           - list of freq, indicator vector and matrices for ch data created by process.ch 
#'# Phi.dm         - list of design matrices; a dm for each capture history 
#'# p.dm           - list of design matrices; a dm for each capture history 
#'# debug          - if TRUE show iterations with par and -2lnl 
#'# time.intervals - intervals of time between occasions 
#'# Value: -LnL using
#'#################################################################################
#'cjs.lnl=function(par,model_data,Phi.links=NULL,p.links=NULL,debug=FALSE,all=FALSE) {
#'	if(debug)cat("\npar = ",par)
#'	#extract Phi and p parameters from par vector 
#'	nphi=ncol(model_data$Phi.dm)
#'	np=ncol(model_data$p.dm) 
#'	beta.phi=par[1:nphi]
#'	beta.p=par[(nphi+1):(nphi+np)] 
#'	#construct parameter matrices (1 row for each capture history and a column 
#'  #for each occasion)
#'	Phis=get.Phi(beta.phi,model_data$Phi.dm,nocc=ncol(model_data$imat$chmat),
#'                           model_data$imat$Fplus)
#'	if(!is.null(model_data$time.intervals)) 
#'	{
#'		exponent=cbind(rep(1,nrow(Phis)),model_data$time.intervals)
#'		Phis=Phis^exponent 
#'	} 
#'	ps=get.p(beta.p,model_data$p.dm,nocc=ncol(model_data$imat$chmat),
#'            model_data$imat$Fplus)
#'	if(debug)cat("\npar = ",par)
#'	# Compute probability of dying in interval from Phis
#'	M=cbind((1-Phis)[,-1],rep(1,nrow(Phis)))
#'	# compute cummulative survival from release across each subsequent time
#'	# and the cummulative probability for detection (capture) across each time
#'	Phi.cumprod=1-model_data$imat$Fplus + Phis*model_data$imat$Fplus 
#'	cump=(1-model_data$imat$Fplus)+model_data$imat$Fplus*
#'          (model_data$imat$chmat*ps+(1-model_data$imat$chmat)*(1-ps))
#'	for (i in 2:ncol(cump))
#'	{
#'		Phi.cumprod[,i]=Phi.cumprod[,i-1]*Phi.cumprod[,i] 
#'		cump[,i]=cump[,i-1]*cump[,i] 
#'	}
#'	# compute prob of capture-history
#'	pch=rowSums(model_data$imat$L*M*Phi.cumprod*cump)
#'	lnl=-sum(model_data$imat$freq*log(pch))
#'	if(debug)cat("\n-2lnl = ",2*lnl) 
#'	return(lnl) 
#'} 
#'}
#' 
#' @param par vector of parameter values
#' @param model_data a list that contains: 1)imat-list of vectors and matrices constructed by
#' \code{\link{process.ch}} from the capture history data, 2)Phi.dm design matrix for Phi constructed by \code{\link{create.dm}},
#' 3)p.dm design matrix for p constructed by \code{\link{create.dm}},
#' 4)Phi.fixed matrix with 3 columns: ch number(i), occasion number(j),
#' fixed value(f) to fix phi(i,j)=f, 5)p.fixed matrix with 3 columns: ch number(i), occasion number(j), and
#' 6) time.intervals intervals of time between occasions if not all 1
#' fixed value(f) to fix p(i,j)=f
#' @param Phi.links vector of links for each parameter
#' @param p.links vector of links for each parameter
#' @param debug if TRUE will printout values of \code{par} and function value
#' @param all if TRUE, returns entire list rather than just lnl; can be used to
#' extract reals
#' @param cjsenv environment for cjs to update iteration counter
#' @return either -log likelihood value if \code{all=FALSE} or the entire
#' list contents of the call to the FORTRAN subroutine if \code{all=TRUE}. The
#' latter is used from \code{\link{cjs}} after optimization to extract the real
#' parameter estimates at the final beta values.
#' @author Jeff Laake
#' @references Pledger, S., K. H. Pollock, et al. (2003). Open
#' capture-recapture models with heterogeneity: I. Cormack-Jolly-Seber model.
#' Biometrics 59(4):786-794.
#' 
# cjs.lnl=function(par,model_data,Phi.links=NULL,p.links=NULL,debug=FALSE,all=FALSE,cjsenv,nc,cl)
cjs.lnl=function(par,model_data,Phi.links=NULL,p.links=NULL,debug=FALSE,all=FALSE,cjsenv)
{
	f_eval=get("markedfunc_eval",envir=cjsenv)+1
	assign("markedfunc_eval", f_eval, envir = cjsenv)
	if(debug)cat("par = ",par,"\n")
	nocc=model_data$imat$nocc
	nphi=ncol(model_data$Phi.dm)
	np=ncol(model_data$p.dm)
	beta.phi=par[1:nphi]
	beta.p=par[(nphi+1):(nphi+np)]
	Phibeta=as.vector(model_data$Phi.dm%*%beta.phi)
#	if(length(Phi.links>0))
#	{
#		Phibeta[Phi.links]=as.vector((sin(Phi.dm[Phi.links,,drop=FALSE]%*%beta.phi)+1)/2)
#	    Phibeta[Phi.links]=log(1/(1/Phibeta[Phi.links]-1))
#	}
	Phibeta=matrix(Phibeta,ncol=nocc-1,nrow=nrow(model_data$Phi.dm)/(nocc-1),byrow=TRUE)
	pbeta=as.vector(model_data$p.dm%*%beta.p)
#	if(length(p.links>0))
#	{
#	   pbeta[p.links]=as.vector((sin(p.dm[p.links,,drop=FALSE]%*%beta.p)+1)/2)
#	   pbeta[p.links]=log(1/(1/pbeta[p.links]-1))
#    }
	pbeta=matrix(pbeta,ncol=nocc-1,nrow=nrow(model_data$p.dm)/(nocc-1),byrow=TRUE)	
#   put Fortran code in a separate function with arguments shown and vector of indices
#   use parLapply with that function and its arguments and index matrix
#   get list - use do.call and sum lnl and use c for p0.
# attempt to parallelize code
#	split_vec=function(nrow,nc)
#	{
#		nx=floor(nrow/nc)
#		return(cbind(seq(1,nrow,nx)[1:nc],
#						c(seq(nx,nrow,nx)[1:(nc-1)],nrow)))
#	}   
#	if(nc>1 & nrow(Phibeta)>5*nc)
#	{
#	   res=parLapply(cl,1:nc, cjs.parallel,chmat=model_data$imat$chmat,Phibeta=Phibeta,pbeta=pbeta,
#			first=model_data$imat$first,last=model_data$imat$last,freq=model_data$imat$freq,
#			loc=model_data$imat$loc,Phi.fixed=model_data$Phi.fixed,p.fixed=model_data$p.fixed,
#			time.intervals=model_data$time.intervals,vec=split_vec(nrow(Phibeta),nc))
#	   lnl=list(lnl=sum(unlist(lapply(res,function(x) x$lnl))))
#	   lnl$p0=unlist(lapply(res,function(x) x$p0))
#   } else
  	   lnl=.Fortran("cjs",as.double(model_data$imat$chmat),as.double(Phibeta),as.double(pbeta),
			as.double(model_data$imat$first),as.double(model_data$imat$last),as.double(model_data$imat$freq),
			as.integer(model_data$imat$loc),as.double(model_data$Phi.fixed),as.double(model_data$p.fixed),
			as.double(model_data$time.intervals),as.integer(nrow(model_data$imat$chmat)),           
			as.integer(ncol(model_data$imat$chmat)),as.integer(nrow(model_data$Phi.fixed)),
			as.integer(nrow(model_data$p.fixed)),lnl=double(1),p0=double(nrow(model_data$imat$chmat)),PACKAGE="marked")
	if(debug)
	{
		cat("-2lnl = ",2*lnl$lnl,"\n")
	} else
	if((f_eval-100*floor(f_eval/100))==0)
	{
	    cat("\r Number of evaluations: ",f_eval," -2lnl:",formatC(2*lnl$lnl,digits=10))
		flush.console()
	}	
	if(all)
		return(lnl)
	else
		return(lnl$lnl)
}
# attempt to parallelize the code
#cjs.parallel=function(i,chmat,Phibeta,pbeta,first,last,freq,loc,Phi.fixed,p.fixed,time.intervals,vec)
#{
#ir=vec[i,1]:vec[i,2]
#lnl=.Fortran("cjs",as.double(chmat[ir,]),as.double(Phibeta[ir,]),as.double(pbeta[ir,]),
#		as.double(first[ir]),as.double(last[ir]),as.double(freq[ir]),
#		as.integer(loc[ir]),as.double(Phi.fixed),as.double(p.fixed),
#		as.double(time.intervals[ir,]),as.integer(nrow(chmat[ir,])),           
#		as.integer(ncol(chmat[ir,])),as.integer(nrow(Phi.fixed)),
#		as.integer(nrow(p.fixed)),lnl=double(1),p0=double(nrow(chmat[ir,])),PACKAGE="marked")
#return(lnl)
#}
