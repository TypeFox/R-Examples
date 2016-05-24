"CovE.n" <-
function(X,y,u,theta,sigma) {
n     <- nrow(X); np <- ncol(X)
rs    <- as.vector(y-X%*%as.matrix(theta))/sigma
wi    <- (abs(rs)<u)*1
sumwi <- sum(wi)
#XtX  <- (t(X)%*%diag(wi)%*%X)/sumwi
#xbar <- apply(wi*X,2,mean)*(n/sumwi)
  wiX <- wi*X 
  XtX <- (t(X) %*% wiX)/sumwi 
 xbar <- apply(wiX, 2, mean) * (n/sumwi)  
U     <- min(u,10); L <- -U
invM0 <- invM2.n(U,theta,sigma,rs,wi,XtX,xbar,estim="SI")$Minv 
invM1 <- invM2.n(U,theta,sigma,rs,wi,XtX,xbar,estim="TMLI")$Minv
xk    <- 1.5477
E.S   <- E.Smat.n(XtX,xbar)
CV0   <- invM0%*%E.S%*%t(invM0)
beta  <- Beta.n(U)
MS0   <- invM0[np+1,,drop=FALSE]
MT0   <- invM0[1:np,,drop=FALSE]

q1    <- (U-L)*dnorm(U)/sigma
#q2   <- 0; q4    <- 0
q5    <- (U^3-L^3)*dnorm(U)/sigma - beta*(U*dnorm(U)-L*dnorm(L))/sigma
E0    <- E0vect.n(xbar)
ku    <- min(xk,U)
E0l   <- E0.vect.n(xbar,ku,U)
E1    <- E1vect.n(xbar,ku,U)
E2    <- E2vect.n(xbar,ku,u)
Q1    <- Q1mat.n(XtX,xbar,E.S,E1,MT0,MS0,U,q1)
Q2    <- Q2vect.n(XtX,xbar,E.S,E0l,E1,E2,MT0,MS0,q1,beta,q5)
Q3    <- Q3sca.n(xbar,E.S,E0l,E2,MT0,MS0,U,beta,q5)
E.Q   <- matrix(0,ncol=np+1,nrow=np+1)
E.Q[1:np,1:np] <- Q1;    E.Q[1:np,np+1] <- Q2
E.Q[np+1,1:np] <- t(Q2); E.Q[np+1,np+1] <- Q3
CV1   <- invM1%*%E.Q%*%t(invM1)
list(CV0=CV0/n, CV1=CV1/n, XtX=XtX, xbar=xbar)}

