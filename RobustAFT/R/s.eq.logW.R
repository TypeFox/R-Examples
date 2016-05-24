S.eq.logW <-
function(X,y,N,delta,sigma,S.ini,b,beta,gama,MAXIT,TOL,mu0=-0.1351788,
             ipsi=4,xk=1.717816,ialg=3) {
X     <- as.matrix(X)
n     <- length(y); p <- ncol(X); mdx   <- nrow(X)
s0    <- 1; cnst <- c(0,0); tint <- 2; meth <- 1
sig5  <- sigma; if (length(sigma)==1) sig5 <- rep(sigma,N)
theta <- rs <- yy <- dd <- sz <- sw <- it <- rep(0,n)
sx    <- matrix(0,nrow=n,ncol=p)
sgama <- sbeta <- rep(0,p)
dfcomn2(ipsi=ipsi,xk=xk)

f.res <- .Fortran("sigama",
x=to.single(X),y=to.single(y),delta=to.single(delta),sig=to.single(S.ini),
mu0=to.single(mu0),s0=to.single(s0),ipsi=to.integer(ipsi),xk=to.single(xk),
b=to.single(b),beta=to.single(beta),gamma=to.single(gama),
cnst=to.single(cnst),n=to.integer(n),np=to.integer(p),ns=to.integer(N),
mdx=to.integer(mdx),lint=to.integer(tint),meth=to.integer(meth),
ialg=to.integer(ialg),maxit=to.integer(MAXIT),tol=to.single(TOL),nit=integer(1),
sigma=single(N),theta=to.single(theta),rs=to.single(rs),yy=to.single(yy),
dd=to.single(dd),sbeta=to.single(sbeta),sgama=to.single(sgama),sx=to.single(sx),
sz=to.single(sz),sw=to.single(sw),sig5=to.single(sig5),it=to.integer(it),
mes2=integer(4))
list(S=f.res$sigma,nit=f.res$nit,mes2=f.res$mes2)}

s.eq.logW <- function(X,y,N,delta,S.ini,b,beta,MAXIT,TOL,mu0= -0.1351788,s0=1,
                 ipsi=4,xk=1.717816,ialg=3,meth=4) {
X     <- as.matrix(X)
n     <- length(y); p <- ncol(X); mdx   <- nrow(X)
cnst  <- c(0,0)
sig5  <- S.ini/s0; if (length(sig5)==1) sig5 <- rep(S.ini/s0,N)
tint  <- 2
theta <- rs <- yy <- dd <- sz <- sw <- it <- rep(0,n)
sx    <- matrix(0,nrow=n,ncol=p)
sgama <- sbeta <- rep(0,p)
gama  <- matrix(0,ncol=p,nrow=N)
dfcomn2(ipsi=ipsi,xk=xk)

f.res <- .Fortran("sigama",
x=to.single(X),y=to.single(y),delta=to.single(delta),sig=to.single(S.ini),
mu0=to.single(mu0),s0=to.single(s0),ipsi=to.integer(ipsi),xk=to.single(xk),
b=to.single(b),beta=to.single(beta),gamma=to.single(gama),
cnst=to.single(cnst),n=to.integer(n),np=to.integer(p),ns=to.integer(N),
mdx=to.integer(mdx),lint=to.integer(tint),meth=to.integer(meth),
ialg=to.integer(ialg),maxit=to.integer(MAXIT),tol=to.single(TOL),nit=integer(1),
sigma=single(N),theta=to.single(theta),rs=to.single(rs),yy=to.single(yy),
dd=to.single(dd),sbeta=to.single(sbeta),sgama=to.single(sgama),sx=to.single(sx),
sz=to.single(sz),sw=to.single(sw),sig5=to.single(sig5),it=to.integer(it),
mes2=integer(4))
list(S=f.res$sigma,nit=f.res$nit,mes2=f.res$mes2)}


