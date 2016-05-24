"CovE.w1" <-
function(y,l,u,theta,sigma) {
n   <- length(y); np <- 1
rs    <- (y-theta)/sigma
wi    <- ((l<rs)&(rs<u))*1
sumwi <- sum(wi)
XtX   <- as.matrix(1)
xbar  <- 1
l     <- max(l,-25); u <- min(u,4)
invM0 <- invM2.w1(l,u,theta,sigma,rs,wi,estim="SI")$Minv 
invM1 <- invM2.w1(l,u,theta,sigma,rs,wi,estim="TMLI")$Minv
xk    <- 1.717817
E.S   <- E.Smat(XtX,xbar)
CV0   <- invM0%*%E.S%*%t(invM0)
beta  <- Beta.w(l,u)
MS0   <- invM0[np+1,,drop=FALSE]
MT0   <- invM0[1:np,,drop=FALSE]
MT0[1,] <- MT0[1,]+0.1352*MS0
q1    <- (exp(u)-exp(l))*ezez(u)/sigma
q2    <- (u*exp(u)-u-l*exp(l)+l)*ezez(u)/sigma
q4    <- (u*(exp(u)-1)-l*(exp(l)-1))*ezez(u)/sigma
q5    <- (u*u*(exp(u)-1)-l*l*(exp(l)-1))*ezez(u)/sigma -
          beta*(u*ezez(u)-l*ezez(l))/sigma
E0    <- E0vect(xbar)
ku    <- min(xk,u); kl <- max(-xk,l)
E0l   <- E0.vect(xbar,kl,ku,l,u)
E1    <- E1vect(xbar,kl,ku,l,u)
E2    <- E2vect(xbar,kl,ku,l,u)
Q1    <- Q1mat(XtX,xbar,E.S,E1,MT0,MS0,l,u,q1,q2)
Q2    <- Q2vect(XtX,xbar,E.S,E0l,E1,E2,MT0,MS0,l,u,q1,q2,beta,q4,q5)
Q3    <- Q3sca(xbar,E.S,E0l,E2,MT0,MS0,l,u,beta,q4,q5) 
E.Q   <- matrix(0,ncol=np+1,nrow=np+1)
E.Q[1:np,1:np] <- Q1;    E.Q[1:np,np+1] <- Q2
E.Q[np+1,1:np] <- t(Q2); E.Q[np+1,np+1] <- Q3
CV1   <- invM1%*%E.Q%*%t(invM1)
list(CV0=CV0/n, CV1=CV1/n)}

