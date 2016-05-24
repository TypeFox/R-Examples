"TML.logWeibull" <-
function(X,y,cu=1.855356,initial=c("S","input"),otp=c("adaptive","fixed"),
                  cov=c("no","parametric","nonparametric"),input=NULL,iv=1,nrep=0,
                  seed=1313,maxit=100,tol=0.0001,gam=0.2,nitmon=FALSE,fastS=TRUE,...){
#
# Main function for TML-estimates; Log-Weibull case
#
X   <- as.matrix(X)         
n   <- length(y); np <- ncol(X); cl <- Izero(cu)
ips <- 2; xk <- 1.717817; beta <- 0.5
namat <- matrix(NA,nrow=np,ncol=np)
nares <- list(th0=NA,v0=NA,th1=rep(NA,np),v1=NA,tl=NA,tu=NA,CV0=namat,V0=NA,CV1=namat,V1=NA,
         alpha=NA,tn=NA,beta=NA,wi=rep(NA,n))
if (all(X[1,]!=1)) {cat("First column of the X matrix must be ones!\n"); return(nares)}
# Step 1: initial high bdp estimate
if (initial=="S")     if (fastS) {
                          set.seed(seed)
                          zctrl <- lmrob.control(...)                         
                          z     <- lmrob.S(X,y,zctrl); th0 <- z$coef; v0 <- z$scale}
                       else { 
                       if (np <= 2 & n <= 500) iopt <- 3 else iopt <- 1; if (nrep!=0) iopt <- 2
#                      z   <- S.hysest(X,y,iopt,nrep,np+1,ips,xk,beta,time=F,seed=seed)
                       z   <- dfcomn2(ipsi=4,xk=xk,beta=beta)
                       z   <- hysest(X,y,nq=np+1,iopt=iopt,intch=1,nrep=nrep,tols=tol,tolr=tol,iseed=seed)
                       th0 <- z$theta[1:np]; v0  <- z$smin
                       b0  <- -0.1352; th0[1] <- th0[1]-b0*v0}
if (initial=="input") {z  <- input; v0 <- z$v; th0 <- z$tau}
nares <- list(th0=th0,v0=v0,th1=rep(NA,np),v1=NA,tl=NA,tu=NA,CV0=namat,V0=NA,CV1=namat,V1=NA,
         alpha=NA,tn=NA,beta=NA,wi=rep(NA,n))
# Step 2: rejection rule
re    <- y-as.vector(X%*%as.matrix(th0)); rs <- re/v0
tp    <- adaptw(sort(rs),cl,cu,otp); if (is.na(tp$tu)) return(nares)
wi    <- tPsiw(rs,tp$tl,tp$tu)
yr    <- y[wi!=0 | rs==0]
Xr    <- X[wi!=0 | rs==0,,drop=FALSE]
tp$tn <- length(yr)
# Step 3: ML-estimate on retained observations
z     <- MLwp(Xr,yr,th0,v0,iv,n,tp,gamm=gam,maxit,tol,nitmon)
res   <- list(th0=th0,v0=v0,th1=z$th1,v1=z$v1,nit=z$nit,tl=tp$tl,tu=tp$tu,
         alpha=tp$alpha,tn=tp$tn,beta=tp$beta,wi=(wi!=0)*1)
if (cov!="no") {l <- cl; u <- cu; if (otp=="adaptive") {l <- tp$tl; u <- tp$tu} 
# if (cov=="halfparametric") K <- Cov2.w(X,y,l,u,z$th1,z$v1,opt="integrals")
  if (cov=="nonparametric") K <- Cov2.w(X,y,l,u,z$th1,z$v1,opt="averages")
  if (cov=="parametric"   ) K <- CovE.w(X,y,l,u,z$th1,z$v1)
  res  <- c(res,list(CV0=K$CV0,CV1=K$CV1))}
res}

