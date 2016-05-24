
np.gcv <- function(data=data, h.seq=NULL, num.h=50, estimator="NW", kernel="quadratic")
{
  
if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data) != 2)  stop("data must have 2 columns: y and t")

if ( (!is.null(h.seq)) && (sum(is.na(h.seq))  != 0) ) stop ("h.seq must be numeric")
if ( (!is.null(h.seq)) && (any(h.seq<=0)) ) stop ("h.seq must contain one ore more positive values")

if (is.null(num.h))   stop ("num.h must not be NULL") 
if (length(num.h) !=1)  stop ("num.h must be an only value")
if (!is.numeric(num.h))   stop ("num.h must be numeric") 
if (num.h<=0)  stop ("num.h must be a positive value") 
  
if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")
  
if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")


kernel.function <- get(kernel)       
n <- nrow(data)
y <- data[, 1]
t <- data[, 2]
    
  
if (is.null(h.seq)) {
  a <- as.matrix(abs(outer(t, t,"-")))
  for (i in 1:n) {a[i,i] <- -1000}
  a <- as.vector(a[a!=-1000])
  
  h.min <- quantile(a,0.05)
  h.max <- (max(t)-min(t))*0.25
  
  h.seq <- seq(h.min, h.max, length.out=num.h) 
} 
else num.h <- length(h.seq)


GCV <- rep(0, num.h)
h.GCV <- data.frame(h=0)

  
W.g <- function(t=t, g=NULL, estimator=estimator, kernel.function=kernel.function)
{
    
  if (estimator=="NW") {
      
    Zmat <- outer(t, t, "-")
    Umat <- Zmat/g
      
    Kmat <- kernel.function(Umat)
    Kmat[(Kmat<0)] <- 0 
      
    S0 <- apply(Kmat, 1, sum)
    Kmat <- Kmat/S0
      
  }      
    
  else if (estimator=="LLP") {
      
    Zmat <- outer(t, t, "-")
    Umat <- Zmat/g
    
    Kmat <- kernel.function(Umat)
    Kmat[(Kmat<0)] <- 0
    
    S0 <- apply(Kmat, 1, sum)
    S1 <- apply(Kmat*Zmat, 1, sum)
    S2 <- apply(Kmat*Zmat^2, 1, sum)
    
    Kmat <- Kmat * (S2 - Zmat*S1)/(S0*S2 - S1^2)
    
  }
}
  
  
for (i in 1:num.h) {
      
  A.h <- W.g(t=t, g=h.seq[i], estimator=estimator, kernel.function=kernel.function)
      
  RSS.h <- sum(((diag(n) - A.h) %*% y)^2)/n
      
  GCV[i] <- RSS.h / (1 - sum(diag(A.h))/n )^2
      
} # for i
    
  

index.GCV <- order(GCV)[1]
  
h.GCV <- c(h.seq[index.GCV])
GCV.opt <- GCV[index.GCV]
  


list(h.opt=h.GCV, GCV.opt=GCV.opt, GCV=GCV, h.seq=h.seq)
  
}
