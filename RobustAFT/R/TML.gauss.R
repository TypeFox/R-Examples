"TML.gauss" <-
function(X,y,cu=2.5,initial=c("S","input"),otp=c("adaptive","fixed"),
             cov=c("no","parametric","nonparametric"),input=NULL,iv=1,nrep=0,
             tol=0.0001,seed=1313,fastS=TRUE,...){
#
# Main function for TML-estimates; Gaussian case
#                
X   <- as.matrix(X)         
n   <- length(y); np <- ncol(X); cl <- -cu
ips <- 2; xk <- 1.54764; beta <- 0.5 
# Step 1: Initial high bdp estimate
if (initial=="S")     {if (fastS) {
                          set.seed(seed)
                          zctrl <- lmrob.control(...)                         
                          z     <- lmrob.S(X,y,zctrl); th0 <- z$coef; v0 <- z$scale}
                       else {                              
                          z  <- dfcomn2(ipsi=4,xk=xk,beta=beta)
                          if (np <= 2 & n <= 500) iopt <- 3 else iopt <- 1; if (nrep!=0) iopt <- 2
                          z  <- hysest(X,y,np+1,iopt=iopt,intch=1,nrep=nrep,tols=tol,tolr=tol,iseed=seed)
                         th0 <- z$theta[1:np]; v0  <- z$smin}
}
if (initial=="input") {z   <- input
                       th0 <- z$lambda;      v0 <- z$sigma}
namat <- matrix(NA,nrow=np,ncol=np)
nares <- list(th0=th0,v0=v0,th1=rep(NA,np),v1=NA,tl=NA,tu=NA,CV0=namat,V0=NA,CV1=namat,V1=NA,
         alpha=NA,tn=NA,beta=NA,wi=rep(NA,n))
# Step 2: rejection rule
re    <- y-as.vector(X%*%as.matrix(th0)); rs <- re/v0
tp    <- adaptn(sort(rs),cl,cu,option=otp); if (is.na(tp$tu)) return(nares)
wi    <- tPsin(rs,tp$tl,tp$tu)
yr    <- y[wi!=0 | rs==0]
Xr    <- X[wi!=0 | rs==0,,drop=FALSE]
tp$tn <- length(yr)
# Step 3: ML-estimate on retained observations
z     <- MLnp(Xr,yr,iv,tp)
res   <- list(th0=th0,v0=v0,th1=z$th1,v1=z$v1,tl=tp$tl,tu=tp$tu,
         alpha=tp$alpha,tn=tp$tn,beta=tp$beta,wi=(wi!=0)*1)
if (cov!="no") {l <- cl; u <- cu; if (otp=="adaptive") {l <- tp$tl; u <- tp$tu} 
# if (cov=="halfparametric") K <- Cov2.n(X,y,u,z$th1,z$v1,opt="integrals",xk=xk)
  if (cov=="nonparametric")  K <- Cov2.n(X,y,u,z$th1,z$v1,opt="averages", xk=xk)
  if (cov=="parametric"   )  K <- CovE.n(X,y,u,z$th1,z$v1)
res  <- c(res,list(CV0=K$CV0,CV1=K$CV1))}
res}

