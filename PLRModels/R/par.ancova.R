
par.ancova <- function(data=data, time.series=FALSE, Var.Cov.eps=NULL, p.max=3, 
                       q.max=3, ic="BIC", num.lb=10, alpha=0.05) 
{

if (!is.array(data))  stop("data must be an array")
if (ncol(data) < 2)  stop("data must have at least 2 columns")

if (!is.logical(time.series)) stop("time.series must be logical")

if ( (!is.null(Var.Cov.eps)) && (sum(is.na(Var.Cov.eps)) != 0) ) stop("Var.Cov.eps must have numeric values")
if ( (!is.null(Var.Cov.eps)) && (!is.array(Var.Cov.eps)) )  stop("Var.Cov.eps must be an array")
if ( (!is.null(Var.Cov.eps)) && ( (ncol(Var.Cov.eps) != nrow(data)) | (nrow(Var.Cov.eps) != nrow(data)) | (dim(Var.Cov.eps)[3] != dim(data)[3]) ) ) stop("Var.Cov.eps must have dimension n x n x L")
if (!is.null(Var.Cov.eps)) {for (k in 1:(dim(Var.Cov.eps)[3]) ) if (any(t(Var.Cov.eps[,,k]) != Var.Cov.eps[,,k])) stop("Var.Cov.eps must be symmetric") }

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
  

  
L <- dim(data)[3]
n <- nrow(data)
p <- ncol(data)-1
  

BETA <- 1:(L*p)  
A <- matrix(0,p*L,p*L)
beta.est <- array(0,c(p,1,L))
B <- matrix(0, p*(L-1), p*L)


B[row(B)==col(B)] <- 1
I <- matrix(0,p,p)
I[row(I)==col(I)] <- 1
for (k in 1:(L-1)) B[(k-1)*p+(1:p),(L-1)*p+(1:p)] <- -I
  


for (k in 1:L)	{
    
   Y <- data[, 1, k]
   X <- data[, -1, k]
   if (!is.matrix(X))  X <- as.matrix(X)
    
   X.X <- t(X)%*%X
   X.Y <- t(X)%*%Y
   beta.est[,1,k] <- symsolve(X.X, X.Y)
        
   BETA[(k-1)*p+(1:p)] <- beta.est[,1,k]
       
   X.X.1 <- solve(X.X)	

   if (!is.null(Var.Cov.eps)) V.eps <- Var.Cov.eps[,,k]

   else {
   		
       eps <- Y - X%*%beta.est[,1,k]
 		
		   if (!time.series) {       
   			  var.eps <- var(eps) 
   			  var.eps <- as.numeric(var.eps)	
  			  V.eps <- diag(var.eps,n,n)
			  }
	
    		else {
          V.eps <- matrix(NA, n, n)
          Var.Cov.mat <- var.cov.matrix(x=eps, n=n, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)

          V.eps <- Var.Cov.mat[[1]]
          pv.Box.test <- Var.Cov.mat[[2]]
          pv.t.test <- Var.Cov.mat[[3]]
          ar.ma <- Var.Cov.mat[[4]]
    	}

   }

   A.k <- n*X.X.1%*%t(X)%*%V.eps%*%X%*%X.X.1
    
   A[(k-1)*p+(1:p) , (k-1)*p+(1:p)] <- A.k
    
} # for k
  

Q.beta  <- n*t(B%*%BETA)%*%solve(B%*%A%*%t(B))%*%(B%*%BETA)
p.v <- 1-pchisq(Q.beta, df=p*(L-1), ncp=0)


if ((is.null(Var.Cov.eps)) && (time.series)) list(par.ancova=data.frame(Q.beta=Q.beta, p.value=p.v), pv.Box.test=pv.Box.test, pv.t.test=pv.t.test, ar.ma=ar.ma)
else list(par.ancova=data.frame(Q.beta=Q.beta, p.value=p.v))

}   

