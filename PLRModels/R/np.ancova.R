
np.ancova <- function(data=data, h.seq=NULL, w= NULL, estimator="NW", kernel="quadratic",
                      time.series=FALSE, Tau.eps=NULL, h0=NULL, lag.max=50, p.max=3, 
                      q.max=3, ic="BIC", num.lb=10, alpha=0.05)
{  
  
if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data)<3)  stop("data must have at least 3 columns")

if ( (!is.null(h.seq)) && (sum(is.na(h.seq))  != 0) ) stop ("h.seq must be numeric")
if ( (!is.null(h.seq)) && (any(h.seq<=0)) ) stop ("h.seq must contain one ore more positive values")
 
if ( (!is.null(w)) && (sum(is.na(w) )  != 0) ) stop ("w must be numeric")
if ( (!is.null(w)) && (!is.vector(w)) ) stop ("w must be a vector")
if ( (!is.null(w)) && (length(w)!=2) ) stop("w must be a vector of length 2")
if ( (!is.null(w)) && (any(w<0)) ) stop ("w must contain two positive values")
if ( (!is.null(w)) && (w[1]>w[2]) ) stop("w[2] must be greater than w[1]")

if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")

if ((kernel=="quadratic") | (kernel=="Epanechnikov") ) {const1 <- 0.6; const2 <- 0.4337}
else if (kernel=="triweight") {const1 <- 0.8159; const2 <- 0.5879}
else if (kernel=="gaussian") {const1 <- 0.2821; const2 <- 0.1995}
else if (kernel=="uniform") {const1 <- 0.5; const2 <- 0.3333}  
else stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")

if (!is.logical(time.series)) stop("time.series must be logical")

if ( (!is.null(Tau.eps)) && (sum(is.na(Tau.eps))  != 0) ) stop("Tau.eps must have numeric values")
if ( (!is.null(Tau.eps)) && (!is.vector(Tau.eps)) )  stop("Tau.eps must be a vector")
if ( (!is.null(Tau.eps)) && (length(Tau.eps) != (ncol(data)-1)) ) stop ("Tau.eps must have length equals to ncol(data)-1")

if ( (!is.null(h0)) && (length(h0) !=1) ) stop ("h0 must be an only value") 
if ( (!is.null(h0)) && (!is.numeric(h0)) )  stop ("h0 must be numeric")
if ( (!is.null(h0)) && (h0<=0) ) stop ("h0 must be a positive value") 

if (is.null(lag.max))   stop ("lag.max must not be NULL") 
if (length(lag.max) !=1)  stop ("lag.max must be an only value")
if (!is.numeric(lag.max))   stop ("lag.max must be numeric") 
if (lag.max<0)  stop ("lag.max must be a positive value") 

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
L <- ncol(data)-1
y <- data[, 1:L]
t <- data[, L+1]



dif <- t[2]-t[1]

for (i in 1:(n-1)) {
  
  dist <- t[i+1]-t[i]
  if (abs(dist - dif) < 2e-15) {}
  else stop("t values must be equidistant")    
  
}
  
  

x0 <- min(t)
x1 <- max(t)

y0 <- (1-0.5)/n
y1 <- (n-0.5)/n

slope <- (y1 - y0)/(x1-x0)
intercept <- y1-slope*x1

t1 <- intercept + slope*t
                  


if (is.null(h.seq)) h.seq <- (seq(0.05, 0.25, length.out=10))/slope
num.h <- length(h.seq)

if (is.null(h0)) h0 <- 0.25/slope

if (is.null(w)) w <- (-intercept + c(0.1, 0.9))/slope



h0.2 <- slope*h0
h.seq.2 <- slope*h.seq
w.2 <- intercept + slope*w



M.est <- array(0,c(n,L,num.h))
  
for (k in 1:L) {
      
   X <- y[,k]
   if (!is.matrix(X))  X <- as.matrix(X)
   m <- np.est(data=cbind(X, t1), newt=t1, h.seq=h.seq.2, estimator=estimator, kernel=kernel)
   M.est[,k,] <- m
      
} # for k
          
  

if (is.null(Tau.eps)) {

		eps.0 <- matrix(0,n,L)

		for (k in 1:L) {
  
   			X <- y[,k]
  			if (!is.matrix(X))  X <- as.matrix(X)
   			m.0 <- np.est(data=cbind(X, t1), newt=t1, h.seq=h0.2, estimator=estimator, kernel=kernel)
   			eps.0[,k] <- y[,k] - m.0

		} # for k

		if (!time.series) Tau.eps.0 <- apply(eps.0,2,var)

		else  {Var.Cov.sum <- var.cov.sum(X=eps.0, lag.max=lag.max, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)

		Tau.eps.0 <- Var.Cov.sum[[1]]
		pv.Box.test <- Var.Cov.sum[[2]]
		pv.t.test <- Var.Cov.sum[[3]]
    ar.ma <- Var.Cov.sum[[4]]
		}
    
    
}

else Tau.eps.0 <- Tau.eps

weights <- 1:n
weights[((weights < n*w.2[1]+0.5) | (weights > n*w.2[2]+0.5))] <- 0
weights[weights != 0] <- 1
Q.m <- matrix(0,num.h,1)
Q.m.normalised <- matrix(0,num.h,1)
p.value <- matrix(0,num.h,1)
  

for (j in 1:num.h) {
  
  for (k in 2:L) for (s in 1:(k-1)) Q.m[j,] <- Q.m[j,] + sum(weights*(M.est[,k,j] - M.est[,s,j])^2)/n
  
  mean.Q.m <-  sum(Tau.eps.0) * const1 *(L-1)*(w.2[2]-w.2[1]) * (n*h.seq.2[j])^(-1)
  sd.Q.m <- (n^2*h.seq.2[j])^(-0.5) * sqrt(2*const2*(sum(Tau.eps.0)^2 + (L^2-2*L)*sum(Tau.eps.0^2)) * (w.2[2]-w.2[1]))
  
  Q.m.normalised[j,] <- (Q.m[j,] - mean.Q.m)/sd.Q.m   

  p.value[j,] <- 1-pnorm(q=Q.m[j,], mean=mean.Q.m, sd=sd.Q.m)
  
}
  

if ((is.null(Tau.eps)) && (time.series)) list(np.ancova=data.frame(h.seq=h.seq, Q.m=Q.m, Q.m.normalised=Q.m.normalised, p.value=p.value), pv.Box.test=pv.Box.test, pv.t.test=pv.t.test, ar.ma=ar.ma)
else list(np.ancova=data.frame(h.seq=h.seq, Q.m=Q.m, Q.m.normalised=Q.m.normalised, p.value=p.value))
  
}
  
