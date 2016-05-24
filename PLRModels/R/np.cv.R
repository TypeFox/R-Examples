
np.cv <- function(data=data, h.seq=NULL, num.h=50, w=NULL, num.ln=1, ln.0=0, step.ln=2, 
                  estimator="NW", kernel="quadratic")
{
  
if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data) != 2)  stop("data must have 2 columns: y and t")

if ( (!is.null(h.seq)) && (sum(is.na(h.seq))  != 0) ) stop ("h.seq must be numeric")
if ( (!is.null(h.seq)) && (any(h.seq<=0)) ) stop ("h.seq must contain one ore more positive values")

if (is.null(num.h))   stop ("num.h must not be NULL") 
if (length(num.h) !=1)  stop ("num.h must be an only value")
if (!is.numeric(num.h))   stop ("num.h must be numeric") 
if (num.h<=0)  stop ("num.h must be a positive value") 
  
if ( (!is.null(w)) && (sum(is.na(w) )  != 0) ) stop ("w must be numeric")
if ( (!is.null(w)) && (!is.vector(w)) ) stop ("w must be a vector")
if ( (!is.null(w)) && (length(w)!=2) ) stop("w must be a vector of length 2")
if ( (!is.null(w)) && (any(w<0)) ) stop ("w must contain two positive values")
if ( (!is.null(w)) && (w[1]>w[2]) ) stop("w[2] must be greater than w[1]")

if (is.null(num.ln))   stop ("num.ln must not be NULL") 
if (length(num.ln) !=1)  stop ("num.ln must be an only value")
if (!is.numeric(num.ln))   stop ("num.ln must be numeric") 
if (num.ln<=0)  stop ("num.ln must be greater or equal than 0") 

if (is.null(ln.0))   stop ("ln.0 must not be NULL") 
if (length(ln.0) !=1)  stop ("ln.0 must be an only value")
if (!is.numeric(ln.0))   stop ("ln.0 must be numeric") 
if (ln.0<0)  stop ("ln.0 must be a positive value") 

if (is.null(step.ln))   stop ("step.ln must not be NULL") 
if (length(step.ln) !=1)  stop ("step.ln must be an only value")
if (!is.numeric(step.ln))   stop ("step.ln must be numeric") 
if (step.ln<=0)  stop ("step.ln must be a positive value") 

if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")
  
if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")



kernel.function <- get(kernel)       
n <- nrow(data)
y <- data[, 1]
t <- data[, 2]
  

if (is.null(w)) w <- c(quantile(t,0.1), quantile(t,0.9))
  

if (is.null(h.seq)) {
  a <- as.matrix(abs(outer(t, t,"-")))
  for (i in 1:n) {a[i,i] <- -1000}
  a <- as.vector(a[a!=-1000])
  
  h.min <- quantile(a,0.05)
  h.max <- (max(t)-min(t))*0.25
  
  h.seq <- seq(h.min, h.max, length.out=num.h) 
}
else num.h <- length(h.seq)


CV <- matrix(0, num.h, num.ln)
CV.opt <- 1:num.ln
h.CV <- data.frame(matrix(0,2,num.ln),row.names=c("ln","h"))
  
ele.seq <- seq(from = ln.0, to =ln.0+(num.ln-1)*step.ln , by = step.ln)

  

if (estimator=="NW") 
    
  for (k in 1:num.ln) {
      
    for (i in 1:n) {
        
      if ((w[1]<=t[i]) & (t[i]<=w[2])) {
          
        diff <- t-t[i]                 
        Zmat <- matrix(rep(diff, num.h), nrow = num.h, byrow = T)
          
        Umat <- Zmat/h.seq
        Kmat <- kernel.function(Umat)
        Kmat[(Kmat<0) | abs(col(Kmat)-i)<=ele.seq[k]] <- 0
          
        S0 <- apply(Kmat, 1, sum)
          
        Kmat <- Kmat/S0
          
        Ymat <- t(y)
        Ymat <- matrix(Ymat, num.h, n, byrow=T)
          
        mhat <- apply(Ymat * Kmat, 1, sum)
          
        CV[, k] <- CV[, k] + (y[i]-mhat)^2
          
      } # if 
    } # i
  } # k

  
else if (estimator=="LLP") 
    
  for (k in 1:num.ln) {
      
    for (i in 1:n)   {
        
      if ((w[1]<=t[i]) & (t[i]<=w[2])) {
          
        diff <- t-t[i]               
        Zmat <- matrix(rep(diff, num.h), nrow = num.h, byrow = T)
          
        Umat <- Zmat/h.seq
        Kmat <- kernel.function(Umat)
        Kmat[(Kmat<0) | abs(col(Kmat)-i)<=ele.seq[k]] <- 0
          
        S0 <- apply(Kmat, 1, sum)
        S1 <- apply(Kmat*Zmat, 1, sum)
        S2 <- apply(Kmat*Zmat^2, 1, sum)
          
        Kmat <- Kmat * (S2 - Zmat*S1)/(S0*S2 - S1^2)
          
        Ymat <- t(y)
        Ymat <- matrix(Ymat, num.h, n, byrow=T)
          
        mhat <- apply(Ymat * Kmat, 1, sum)
                
        CV[, k] <- CV[, k] +(y[i]-mhat)^2
          
      } # if
    } # i
  } # k


  
CV<-CV/n
  
for (k in 1:num.ln) {
    
  index.CV <- order(CV[, k])[1]
    
  h.CV[,k] <- c(ele.seq[k], h.seq[index.CV]) 
  CV.opt[k] <- CV[index.CV, k]
  
}  
  
list(h.opt=h.CV, CV.opt=CV.opt, CV=CV, w=w, h.seq=h.seq)
  
}

