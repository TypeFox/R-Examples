#
# -----------------------------------------------------------------------------
#
rhierBinLogit=
function(Data,Prior,Mcmc){
#
# revision history: 
#	changed 5/12/05 by Rossi to add error checking
#       1/07 removed init.rmultiregfp
#       3/07 added classes
#
# purpose: run binary heterogeneous logit model 
#
# Arguments:
#   Data contains a list of (lgtdata[[i]],Z)
#      lgtdata[[i]]=list(y,X)
#         y is index of brand chosen, y=1 is exp[X'beta]/(1+exp[X'beta])
#         X is a matrix that is n_i x by nvar
#      Z is a matrix of demographic variables nlgt*nz that have been 
#	  mean centered so that the intercept is interpretable
#   Prior contains a list of (nu,V,Deltabar,ADelta)
#      beta_i ~ N(Z%*%Delta,Vbeta)
#      vec(Delta) ~ N(vec(Deltabar),Vbeta (x) ADelta^-1)
#      Vbeta ~ IW(nu,V)
#   Mcmc is a list of (sbeta,R,keep)
#      sbeta is scale factor for RW increment for beta_is
#      R is number of draws
#      keep every keepth draw
#
# Output:
#      a list of Deltadraw (R/keep x nvar x nz), Vbetadraw (R/keep x nvar**2), 
#         llike (R/keep), betadraw is a nlgt x nvar x nz x R/keep array of draws of betas
#         nunits=length(lgtdata)
#
#  define functions needed
#
# ------------------------------------------------------------------------
#
loglike=
function(y,X,beta) {
# function computer log likelihood of data for binomial logit model
# Pr(y=1) = 1 - Pr(y=0) = exp[X'beta]/(1+exp[X'beta])
prob = exp(X%*%beta)/(1+exp(X%*%beta))
prob = prob*y + (1-prob)*(1-y)
sum(log(prob))
}
#
#
#  check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of m,lgtdata, and (possibly) Z")}
  if(is.null(Data$lgtdata)) {pandterm("Requires Data element lgtdata (list of data for each unit)")}
  lgtdata=Data$lgtdata
  nlgt=length(lgtdata)
if(is.null(Data$Z)) { cat("Z not specified -- putting in iota",fill=TRUE); fsh() ; Z=matrix(rep(1,nlgt),ncol=1)}
  else {if (nrow(Data$Z) != nlgt) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number logits ",nlgt))}
      else {Z=Data$Z}}
  nz=ncol(Z)
#
# check lgtdata for validity
#
m=2  # set two choice alternatives for Greg's code
ypooled=NULL
Xpooled=NULL
if(!is.null(lgtdata[[1]]$X)) {oldncol=ncol(lgtdata[[1]]$X)}
for (i in 1:nlgt) 
{
    if(is.null(lgtdata[[i]]$y)) {pandterm(paste("Requires element y of lgtdata[[",i,"]]"))}
    if(is.null(lgtdata[[i]]$X)) {pandterm(paste("Requires element X of lgtdata[[",i,"]]"))}
    ypooled=c(ypooled,lgtdata[[i]]$y)
    nrowX=nrow(lgtdata[[i]]$X)
    if((nrowX) !=length(lgtdata[[i]]$y)) {pandterm(paste("nrow(X) ne length(yi); exception at unit",i))}
    newncol=ncol(lgtdata[[i]]$X)
    if(newncol != oldncol) {pandterm(paste("All X elements must have same # of cols; exception at unit",i))}
    Xpooled=rbind(Xpooled,lgtdata[[i]]$X)
    oldncol=newncol
}
nvar=ncol(Xpooled)
levely=as.numeric(levels(as.factor(ypooled)))
if(length(levely) != m) {pandterm(paste("y takes on ",length(levely)," values -- must be = m"))}
bady=FALSE
for (i in 0:1 )
{
    if(levely[i+1] != i) bady=TRUE
}
cat("Table of Y values pooled over all units",fill=TRUE)
print(table(ypooled))
if (bady) 
  {pandterm("Invalid Y")}
#
# check on prior
#
if(missing(Prior)){
    nu=nvar+3
    V=nu*diag(nvar)
    Deltabar=matrix(rep(0,nz*nvar),ncol=nvar)
    ADelta=.01*diag(nz) }
else {
    if(is.null(Prior$nu)) {nu=nvar+3}  else {nu=Prior$nu}
        if(nu < 1) {pandterm("invalid nu value")}
    if(is.null(Prior$V)) {V=nu*diag(rep(1,nvar))} else {V=Prior$V}
    if(sum(dim(V)==c(nvar,nvar)) !=2) pandterm("Invalid V in prior")
    if(is.null(Prior$ADelta) ) {ADelta=.01*diag(nz)} else {ADelta=Prior$ADelta}
    if(ncol(ADelta) != nz | nrow(ADelta) != nz) {pandterm("ADelta must be nz x nz")}
    if(is.null(Prior$Deltabar) ) {Deltabar=matrix(rep(0,nz*nvar),ncol=nvar)} else {Deltabar=Prior$Deltabar}
}
#
# check on Mcmc
#
if(missing(Mcmc)) 
  {pandterm("Requires Mcmc list argument")}
else 
   { 
    if(is.null(Mcmc$sbeta)) {sbeta=.2} else {sbeta=Mcmc$sbeta}
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    }
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Attempting MCMC Inference for Hierarchical Binary Logit:",fill=TRUE)
cat(paste("  ",nvar," variables in X"),fill=TRUE)
cat(paste("  ",nz," variables in Z"),fill=TRUE)
cat(paste("   for ",nlgt," cross-sectional units"),fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("nu =",nu,fill=TRUE)
cat("V ",fill=TRUE)
print(V)
cat("Deltabar",fill=TRUE)
print(Deltabar)
cat("ADelta",fill=TRUE)
print(ADelta)
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat(paste("sbeta=",round(sbeta,3)," R= ",R," keep= ",keep),fill=TRUE)
cat("",fill=TRUE)

nlgt=length(lgtdata)
nvar=ncol(lgtdata[[1]]$X)
nz=ncol(Z)



#
# initialize storage for draws
#
Vbetadraw=matrix(double(floor(R/keep)*nvar*nvar),ncol=nvar*nvar)
betadraw=array(double(floor(R/keep)*nlgt*nvar),dim=c(nlgt,nvar,floor(R/keep)))
Deltadraw=matrix(double(floor(R/keep)*nvar*nz),ncol=nvar*nz)
oldbetas=matrix(double(nlgt*nvar),ncol=nvar)
oldVbeta=diag(nvar)
oldVbetai=diag(nvar)
oldDelta=matrix(double(nvar*nz),ncol=nvar)

betad = array(0,dim=c(nvar))
betan = array(0,dim=c(nvar))
reject = array(0,dim=c(R/keep))
llike=array(0,dim=c(R/keep))

itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min)",fill=TRUE)
fsh()
for (j in 1:R) {
	rej = 0
	logl = 0
	sV = sbeta*oldVbeta
	root=t(chol(sV))

#	Draw B-h|B-bar, V

	for (i in 1:nlgt) {

		betad = oldbetas[i,]
		betan = betad + root%*%rnorm(nvar)
# data		
		lognew = loglike(lgtdata[[i]]$y,lgtdata[[i]]$X,betan)
		logold = loglike(lgtdata[[i]]$y,lgtdata[[i]]$X,betad) 
# heterogeneity
logknew = -.5*(t(betan)-Z[i,]%*%oldDelta) %*% oldVbetai %*% (betan-t(Z[i,]%*%oldDelta))
logkold = -.5*(t(betad)-Z[i,]%*%oldDelta) %*% oldVbetai %*% (betad-t(Z[i,]%*%oldDelta))
# MH step
		alpha = exp(lognew + logknew - logold - logkold)
		if(alpha=="NaN") alpha=-1
		u = runif(n=1,min=0, max=1)
		if(u < alpha) { 
			oldbetas[i,] = betan
			logl = logl + lognew } else {
		 	logl = logl + logold
			rej = rej+1  }
		}
#	Draw B-bar and V as a multivariate regression
	out=rmultireg(oldbetas,Z,Deltabar,ADelta,nu,V)
	oldDelta=out$B
	oldVbeta=out$Sigma
	oldVbetai=chol2inv(chol(oldVbeta))

	if((j%%100)==0) 
          {
           ctime=proc.time()[3]
           timetoend=((ctime-itime)/j)*(R-j)
           cat(" ",j," (",round(timetoend/60,1),")",fill=TRUE)
           fsh() }
	mkeep=j/keep
	if(mkeep*keep == (floor(mkeep)*keep))
          {Deltadraw[mkeep,]=as.vector(oldDelta)
           Vbetadraw[mkeep,]=as.vector(oldVbeta)
           betadraw[,,mkeep]=oldbetas
           llike[mkeep]=logl
           reject[mkeep]=rej/nlgt
          }
}
ctime=proc.time()[3]
cat(" Total Time Elapsed: ",round((ctime-itime)/60,2),fill=TRUE)


attributes(betadraw)$class=c("bayesm.hcoef")
attributes(Deltadraw)$class=c("bayesm.mat","mcmc")
attributes(Deltadraw)$mcpar=c(1,R,keep)
attributes(Vbetadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(Vbetadraw)$mcpar=c(1,R,keep)

return(list(betadraw=betadraw,Vbetadraw=Vbetadraw,Deltadraw=Deltadraw,llike=llike,reject=reject))
}


