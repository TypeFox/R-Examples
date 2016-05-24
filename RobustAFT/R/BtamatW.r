BtamatW <-
function(X,y,delta,N,q,MAXIT,TOL,seed=153)     {            
# Function for determining the matrix of betas; LogWeibull case
n  <- length(y); p <- ncol(X); if (q < p) q <- p
set.seed(seed)
indu  <- (1:n)[delta==1]
inds  <- apply(matrix(rep(indu,N),nrow=N,byrow=TRUE),1,sample,size = q)
intcp <- any(X[,1,drop=TRUE]!= 1)
if (intcp) X <- cbind(1,X)
beta  <- apply(inds, 2, CandidateW, X, y, delta,MAXIT,TOL)
if (p==1) beta <- matrix(beta,ncol=1,nrow=N) else beta <- t(beta)
list(beta = beta)}

