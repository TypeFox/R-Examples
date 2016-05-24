lmmlasso <- function(x,...)
  UseMethod("lmmlasso")

lmmlasso.default <- function(x,y,z=x,grp,weights=NULL,coefInit=NULL,lambda,startValue=1,nonpen=1:dim(z)[[2]],pdMat=c("pdIdent","pdDiag","pdSym"),
                       method="ML",CovOpt=c("nlminb","optimize"),stopSat=TRUE,standardize=TRUE,control=lmmlassoControl(),
                        ranInd=1:dim(z)[[2]],...)
{

# --- Introductory checks and reformulations ---
# ----------------------------------------------

CovOpt <- match.arg(CovOpt)
pdMat <- match.arg(pdMat)
  
# do some checks
if (!is.matrix(x)) stop("x has to be a matrix")
if (any(is.na(x))) stop("Missing values in x not allowed")
if (any(is.na(y))) stop("Missing values in y not allowed")
if (!is.matrix(z)) stop("z has to be a matrix")
if (any(is.na(z))) stop("Missing values in z not allowed")
if (!is.numeric(y)) stop("y has to be of type 'numeric'")
if (nrow(x)!=length(y)) stop("x and y have not correct dimensions")
if (any(x[,1]!=rep(1,dim(x)[[1]]))) stop("first column is not the intercept")
if (length(levels(grp))==1) stop("Only one group. No covariance parameters!")
if (!all(nonpen%in%1:dim(x)[[2]])) stop("Error with the argument nonpen")
if (length(lambda)!=1) stop("lambda can only be one number")
if (lambda<0) stop("regularization parameter must be positive")

# calculate weights (NA: not penalized, 0: drop)
if (missing(weights))
 {
   weights <- rep(1,dim(x)[[2]])
 } else
    {
      if (!length(weights)==dim(x)[[2]]) stop("Weights vector has not length p")

      nonpen <- which(is.na(weights))
      if (any(weights[-nonpen]<0)) stop("Weights must be positive")
      remove <- weights==0
      remove[nonpen] <- FALSE
      
      x <- x[,!remove,drop=FALSE]
      
      rem <- which(weights==0)
      z <- z[,!(ranInd%in%rem),drop=FALSE]
      
      weights <- weights[!remove]
      nonpen <- which(is.na(weights))
      
      ran1 <- logical(dim(x)[[2]])
      ran1[ranInd] <- TRUE
      ran2 <- ran1[!remove]
      ranInd <- which(ran2)
   }

if (standardize)
  {
    xOr <- x
    meanx <- apply(x[,-1],2,mean)
    sdx <- apply(x[,-1],2,sd)
    x <- cbind(1,scale(x[,-1],center=meanx,scale=sdx))
  }

# crucial allocations
grp <- factor(grp)
N <- length(levels(grp)) # N is the number of groups
p <- dim(x)[[2]]         # p is the number of covariates
q <- dim(z)[[2]]         # q is the number of random effects variables
ntot <- length(y)        # ntot is the total number of observations
Q <- q*(q+1)/2           # maximum number of variance components parameters

# save the grouped information as components of a list
yGrp <- split(y,grp)
xGrp <- split.data.frame(x,grp)
zGrp <- split.data.frame(z,grp)
zIdGrp <- mapply(ZIdentity,zGrp)

# --- Calculating or adopting the starting values ---
# ---------------------------------------------------

if (missing(coefInit))
 {
   # --- starting value for beta ---
   set.seed(control$seed)
   if (startValue==1) # 10-fold cv Lasso
     {
       init <- optL1(y,x[,-1],model="linear",fold=10,trace=FALSE)
       betaStart <- c(init$fullfit@unpenalized,init$fullfit@penalized)
       
     } else if (startValue==2) # 10-fold cv Ridge Regression
     {
       init <- optL2(y,x[,-1],model="linear",fold=10,trace=FALSE)
       betaStart <- c(init$fullfit@unpenalized,init$fullfit@penalized)
     } else # Null start
     {
       betaStart <- c(mean(y),rep(0,p-1))
     }
   
   # --- starting values for the variance components parameters ---
   covStart <- covStartingValues(xGrp,yGrp,zGrp,zIdGrp,betaStart,ntot,N)
   sStart <- covStart$sigma
   tStart <- covStart$tau

   if (pdMat=="pdSym")    parsStart <- vecli(tStart*diag(q))
   if (pdMat=="pdDiag")   parsStart <- rep(tStart,q)
   if (pdMat=="pdIdent")  parsStart <- tStart
 
 } else
 {
   betaStart <- coefInit[[1]]
   parsStart <- coefInit[[2]]
   sStart <- coefInit[[3]]
 }
   
if (pdMat=="pdSym")    PsiStart <- crossprod(triang(parsStart,q))
if (pdMat=="pdDiag")   PsiStart <- diag(parsStart^2, nrow = length(parsStart))
if (pdMat=="pdIdent")  PsiStart <- parsStart^2*diag(q)

# --- Calculate objective function for the starting values ---
# ------------------------------------------------------------

lambdaInvGrp <- mapply(LambdaInv,Z=zGrp,ZId=zIdGrp,MoreArgs=list(Psi=PsiStart,sigma=sStart))
ll1 <- 1/2*ntot*log(2*pi)
fctStart <- ObjFunction(xGroup=xGrp,yGroup=yGrp,LGroup=lambdaInvGrp,b=betaStart,weights=weights,lambda=lambda,nonpen=nonpen,ll1=ll1)
if (control$trace>3) cat("fctStart:",fctStart,"\n")


# --- Coordinate Gradient Descent-iteration ---
# ---------------------------------------------

# some necessary allocations:
betaIter <- betaStart
sIter <- sStart
parsIter <- parsStart
PsiIter <- PsiStart
convPar <- crossprod(betaIter)
convCov <- crossprod(c(sStart,parsStart))
fctIter <- convFct <- fctStart
hessian0 <- rep(0,p)
mat0 <- matrix(0,ncol=p,nrow=N)

stopped <- FALSE
doAll <- FALSE
converged <- 0
counterIn <- 0
counter <- 0       # counts the number of outer iterations

while ((control$maxIter > counter) &  (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll ))
 {
  if (control$maxIter==counter+1) {if (control$trace>2) cat("maxIter reached","\n") ; converged <- converged + 1}

  # arguments that change at each iteration
  counter <- counter + 1 ; if (control$trace==2) cat(counter,"...","\n")
  betaIterOld <- betaIter
  fctIterOld <- fctIter
  covIterOld <- c(sIter,parsIter)
  
  # --- optimization w.r.t the fixed effects vector beta ---
  # --------------------------------------------------------

  activeSet <- which(betaIter!=0)
  if ((length(activeSet)>min(p,ntot))&(lambda>0)&(stopSat)&(counter>2)) {stopped <- TRUE ; break}
  
  if (counterIn==0 | counterIn>control$number)
    {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1    
    } else
    {
      doAll <- FALSE
      counterIn <- counterIn+1
    }

  # calculate the hessian matrices for j in the activeSet
  HessIter <- HessIterTrunc <- HessianMatrix(xGroup=xGrp,LGroup=lambdaInvGrp,activeSet=activeSet,N=N,hessian=hessian0,mat=mat0[,activeSet,drop=FALSE])
  HessIter[activeSet] <- pmin(pmax(HessIter[activeSet],control$lower),control$upper)

  LxGrp <- as1(xGrp,lambdaInvGrp,activeSet,N=N)

  ll2 <- nlogdet(LGroup=lambdaInvGrp)
  
  for (j in activeSet)
   {
    cut1 <- as2(x=x,y=y,b=betaIter,j=j,activeSet=activeSet,group=grp,sGroup=LxGrp)
    JinNonpen <- j%in%nonpen

    # optimum can be calculated analytically
    if (HessIterTrunc[j]==HessIter[j])
      {
         if (JinNonpen) {betaIter[j] <- cut1/HessIter[j]} else {betaIter[j] <- SoftThreshold(cut1,lambda/weights[j])/HessIter[j]}
      } else
    
    # optimimum is determined by the armijo rule
      {
        armijo <- ArmijoRule(xGroup=xGrp,yGroup=yGrp,LGroup=lambdaInvGrp,b=betaIter,j=j,cut=cut1,HkOldJ=HessIterTrunc[j],
                             HkJ=HessIter[j],JinNonpen=JinNonpen,lambda=lambda,weights=weights,nonpen=nonpen,
                             ll1=ll1,ll2=ll2,converged=converged,control=control)
        betaIter <- armijo$b
        converged <- armijo$converged
        fctIter <- armijo$fct
      }
    if (control$trace>3) {fctIter <- ObjFunction(xGroup=xGrp,yGroup=yGrp,LGroup=lambdaInvGrp,b=betaIter,weights=weights,
                        lambda=lambda,nonpen=nonpen,ll1=ll1,ll2=ll2) ; print(fctIter) }
   } 

  if (control$trace>3) cat("-----------","\n")

  # --- optimization w.r.t the variance components parameters ---
  # -------------------------------------------------------------
  
  # calculations before the covariance optimization
  activeset <- which(betaIter!=0)
  resGrp <- ResAsSplit(x=x,y=y,b=betaIter,f=grp,activeset=activeset)
  ll4 <- lambda*sum(abs(betaIter[-nonpen])/weights[-nonpen])

  # optimization of the free parameters in \Psi
  if (pdMat=="pdSym")
   { 
     covParOpt <- MLpdSym(zGroup=zGrp,zIdGroup=zIdGrp,resGroup=resGrp,sigma=sIter,pars=parsIter,q=q,
                      thres=control$thres,ll1=ll1,ll4=ll4,trace=control$trace,CovOpt=CovOpt,VarInt=control$VarInt,CovInt=control$CovInt)
     PsiIter <- crossprod(triang(covParOpt$pars,q))
   }
  else if (pdMat=="pdDiag")
    {
     covParOpt <- MLpdDiag(zGroup=zGrp,zIdGroup=zIdGrp,resGroup=resGrp,sigma=sIter,pars=parsIter,q=q,
                       thres=control$thres,ll1=ll1,ll4=ll4,trace=control$trace,CovOpt=CovOpt,VarInt=control$VarInt)
     PsiIter <- diag(covParOpt$pars^2,nrow=length(covParOpt$pars))
    }
  else if (pdMat=="pdIdent")
   {
     covParOpt <- MLpdIdent(zGroup=zGrp,zIdGroup=zIdGrp,resGroup=resGrp,sigma=sIter,pars=parsIter,q=q,
                            thres=control$thres,ll1=ll1,ll4=ll4,trace=control$trace,CovOpt=CovOpt,VarInt=control$VarInt)
     PsiIter <- covParOpt$par^2*diag(q)    
   }
  
  parsIter <- covParOpt$pars
  fctIter <- covParOpt$fct
  if (control$trace>3)   cat("-----------","\n")

  # optimization of the error variance \sigma^2
  covParOpt <- MLsigma(zGroup=zGrp,zIdGroup=zIdGrp,resGroup=resGrp,q=q,ll1=ll1,ll4=ll4,true.sigma=sIter,Psi=PsiIter,
                     trace=control$trace,CovOpt=CovOpt,VarInt=control$VarInt)
  sIter <- covParOpt$sigma
  fctIter <- covParOpt$fct
  
  if (control$trace>3)  cat("-----------","\n") #cat("obj.fct III:",fctIter,"\n")

  lambdaInvGrp <- mapply(LambdaInv,Z=zGrp,ZId=zIdGrp,MoreArgs=list(Psi=PsiIter,sigma=sIter))
  covIter <- c(sIter,parsIter)

  # --- check convergence ---
  convPar <- sqrt(crossprod(betaIter-betaIterOld))/(1+sqrt(crossprod(betaIter)))
  convFct <- abs((fctIterOld-fctIter)/(1+abs(fctIter)))
  convCov <- sqrt(crossprod(covIter-covIterOld))/(1+sqrt(crossprod(covIter)))
  
  if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
  
 }

if (standardize)
  {
    betaIter[-1] <- betaIter[-1]/sdx
    betaIter[1] <- betaIter[1] - sum(meanx*betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr,grp)
  }


# --- prediction of the random effects ---
# ----------------------------------------

biGroup <- u <- list() ; length(biGroup) <- length(u) <- N
Psi <- PsiIter
corPsi <- cov2cor(Psi)

if (pdMat=="pdSym") cholPsi <- t(triang(parsIter,q))
else if (pdMat=="pdDiag") cholPsi <- diag(parsIter,nrow=length(parsIter))
else if (pdMat=="pdIdent") cholPsi <- parsIter*diag(q)

for (i in 1:N)
{
 u[[i]] <- sIter*solve(t(zGrp[[i]]%*%cholPsi)%*%zGrp[[i]]%*%cholPsi+sIter^2*diag(q))%*%t(zGrp[[i]]%*%cholPsi)%*%resGrp[[i]]
 biGroup[[i]] <- 1/sIter*cholPsi%*%u[[i]]
}


#  --- final calculations ---
# ---------------------------

# fitted values and residuals
residGrp <- fittedGrp <- list() ; length(residGrp) <- length(fittedGrp) <- N
for (i in 1:N)
 {
  fittedGrp[[i]] <- xGrp[[i]][,activeSet,drop=FALSE]%*%betaIter[activeSet,drop=FALSE] + zGrp[[i]]%*%biGroup[[i]]
  residGrp[[i]] <- yGrp[[i]]-fittedGrp[[i]]
 }
residuals <- unlist(residGrp)
fitted <- unlist(fittedGrp)

# random effects, sorted per subject
u <- unlist(u) # corresponds to lmer@u
bi <- unlist(biGroup) # unsorted random effects

# random effects, sorted per effect
ranef <- bi[order(rep(1:q,N))] # corresponds to lmer@ranef

# fixed effects without names
fixef <- betaIter
names(fixef) <- NULL

# --- summary information ---
# ---------------------------
npar <- sum(betaIter!=0) + length(c(sIter,parsIter))
logLik <- MLloglik(xGroup=xGrp,yGroup=yGrp,LGroup=lambdaInvGrp,b=betaIter,ntot=ntot,N=N,activeSet=which(betaIter!=0))
deviance <- -2*logLik
aic <- -2* logLik + 2*npar
bic <- -2* logLik + log(ntot)*npar

if (any(parsIter==0)) cat("Redundant covariance parameters.","\n")
if (converged>0) cat("Algorithm does not properly converge.","\n")
if (stopped) {cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.","\n")
                    ; sIter <- parsIter <- nlogLik <- aic <- bic <- NA ; betaIter <- rep(NA,p) ; bi <- fitted <- residuals <- NULL}

out <- list(data=list(x=x,y=y,z=z,grp=grp),weights=weights,coefInit=list(betaStart=betaStart,parsStart=parsStart,sStart=sStart),
                 lambda=lambda,sigma=abs(sIter),pars=parsIter,coefficients=betaIter,random=bi,u=u,ranef=ranef,fixef=fixef,fitted.values=fitted,
                 residuals=residuals,Psi=Psi,corPsi=corPsi,converged=converged,logLik=logLik,npar=npar,deviance=deviance,
                 aic=aic,bic=bic,nonpen=nonpen,counter=counter,pdMat=pdMat,method=method,CovOpt=CovOpt,control=control,call=match.call(),
                 stopped=stopped,ranInd=ranInd,objective=fctIter)

out
structure(out,class="lmmlasso")

}
