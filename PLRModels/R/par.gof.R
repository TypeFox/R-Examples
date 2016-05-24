
par.gof <- function(data=data, beta0=NULL, time.series=FALSE, Var.Cov.eps=NULL, 
                    p.max=3, q.max=3, ic="BIC", num.lb=10, alpha=0.05) 
{
  
if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data) < 2)  stop("data must have at least 2 columns")

if ( (!is.null(beta0)) && (sum(is.na(beta0))  != 0) ) stop("beta0 must have numeric values")
if ( (!is.null(beta0)) && (!is.vector(beta0)) ) stop("beta0 must be a vector")
if ( (!is.null(beta0)) && (length(beta0) != (ncol(data)-1)) ) stop("beta0 must have length equals to ncol(data)-2")

if (!is.logical(time.series)) stop("time.series must be logical")

if ( (!is.null(Var.Cov.eps)) && (sum(is.na(Var.Cov.eps))  != 0) ) stop("Var.Cov.eps must have numeric values")
if ( (!is.null(Var.Cov.eps)) && (!is.matrix(Var.Cov.eps)) )  stop("Var.Cov.eps must be a matrix")
if ( (!is.null(Var.Cov.eps)) && ( (ncol(Var.Cov.eps) != nrow(data)) | (nrow(Var.Cov.eps) != nrow(data)) ) ) stop("Var.Cov.eps must have dimension n x n")
if ( (!is.null(Var.Cov.eps)) && (any(t(Var.Cov.eps) != Var.Cov.eps)  ) ) stop("Var.Cov.eps must be symmetric")

if (is.null(p.max))   stop ("p.max must not be NULL") 
if (length(p.max) !=1)  stop ("p.max must be an only value")
if (!is.numeric(p.max))   stop ("p.max must be numeric") 
if (p.max<0)  stop ("p.max must be a positive value") 

if (is.null(q.max))   stop ("q.max must not be NULL") 
if (length(q.max) !=1)  stop ("q.max must be an only value")
if (!is.numeric(q.max))   stop ("q.max must be numeric") 
if (q.max<0)  stop ("q.max must be a positive value") 

if ( (ic != "BIC") & (ic != "AIC") & (ic != "AICC") )  stop("ic=BIC or ic=AIC or ic=AICC is required")

if (is.null(num.lb))   stop ("num.lb must not be NULL") 
if (length(num.lb) !=1)  stop ("num.lb must be an only value")
if (!is.numeric(num.lb))   stop ("num.lb must be numeric") 
if (num.lb<=0)  stop ("num.lb must be a positive value") 

if (is.null(alpha))   stop ("alpha must not be NULL") 
if (length(alpha) !=1)  stop ("alpha must be an only value")
if (!is.numeric(alpha))   stop ("alpha must be numeric") 
if ( (alpha<0) | (alpha>1) )  stop ("alpha must be between 0 and 1") 



n <- nrow(data)
p <- ncol(data)-1

if (is.null(beta0)) beta0 <- rep(0,length.out=p)
if (!is.matrix(beta0))  beta0 <- as.matrix(beta0)

Y <- data[, 1]
X <- data[, -1]
if (!is.matrix(X))  X <- as.matrix(X) 

X.X <- t(X)%*%X
X.Y <- t(X)%*%Y
beta.est <- symsolve(X.X, X.Y)
        
X.X.1 <- solve(X.X)


if (is.null(Var.Cov.eps)) {
  v.c.eps <- TRUE

	eps <- Y - X%*%beta.est

	if (!time.series) {       
   			var.eps <- var(eps) 
   			var.eps <- as.numeric(var.eps)	
  			Var.Cov.eps <- diag(var.eps,n,n)
	}
	
  else {
	      Var.Cov.eps <- matrix(NA, n, n)
	      Var.Cov.mat <- var.cov.matrix(x=eps, n=n, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)
	
	      Var.Cov.eps <- Var.Cov.mat[[1]]
	      pv.Box.test <- Var.Cov.mat[[2]]
	      pv.t.test <- Var.Cov.mat[[3]]
        ar.ma <- Var.Cov.mat[[4]]
  }
}
else v.c.eps <- FALSE


A <- n*X.X.1%*%t(X)%*%Var.Cov.eps%*%X%*%X.X.1
    
dif <- beta.est - beta0  
Q.beta  <- n*(t(dif) %*% solve(A) %*% dif)
p.v <- 1-pchisq(Q.beta ,df=p,ncp=0)
  

if ((v.c.eps) && (time.series))  list(par.gof=data.frame(Q.beta=Q.beta, p.value=p.v), pv.Box.test=pv.Box.test, pv.t.test=pv.t.test, ar.ma=ar.ma)
else list(par.gof=data.frame(Q.beta=Q.beta, p.value=p.v))
}   

