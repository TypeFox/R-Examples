
plrm.cv <- function(data=data, b.equal.h=TRUE, b.seq=NULL, h.seq=NULL, num.b=NULL, num.h=NULL,
                    w=NULL, num.ln=1, ln.0=0, step.ln=2, estimator="NW", kernel="quadratic")
{

if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data)<3)  stop("data must have at least 3 columns: y, x, t")
  
if (!is.logical(b.equal.h)) stop("b.equal.h must be logical")

if ( (!is.null(b.seq)) && (sum(is.na(b.seq))  != 0) ) stop ("b.seq must be numeric")
if ( (!is.null(b.seq)) && (any(b.seq<=0)) ) stop ("b.seq must contain one ore more positive values")

if ( (!is.null(h.seq)) && (sum(is.na(h.seq))  != 0) ) stop ("h.seq must be numeric")
if ( (!is.null(h.seq)) && (any(h.seq<=0)) ) stop ("h.seq must contain one ore more positive values")

if ( (!is.null(num.b)) && (length(num.b) !=1) ) stop ("num.b must be an only value") 
if ( (!is.null(num.b)) && (!is.numeric(num.b)) )  stop ("num.b must be numeric") 
if ( (!is.null(num.b)) && (num.b<=0) ) stop ("num.b must be a positive value") 

if ( (!is.null(num.h)) && (length(num.h) !=1) ) stop ("num.h must be an only value") 
if ( (!is.null(num.h)) && (!is.numeric(num.h)) )  stop ("num.h must be numeric") 
if ( (!is.null(num.h)) && (num.h<=0) ) stop ("num.h must be a positive value") 

if ( (!is.null(w)) && (sum(is.na(w) )  != 0) ) stop ("w must be numeric")
if ( (!is.null(w)) && (!is.vector(w)) ) stop ("w must be a vector")
if ( (!is.null(w)) && (length(w)!=2) ) stop("w must be a vector of length 2")
if ( (!is.null(w)) && (any(w<0)) ) stop ("w must contain two positive values")
if ( (!is.null(w)) && (w[1]>w[2]) ) stop("w[2] must be greater than w[1]")

if (is.null(num.ln))   stop ("num.ln must not be NULL") 
if (length(num.ln) !=1)  stop ("num.ln must be an only value")
if (!is.numeric(num.ln))   stop ("num.ln must be numeric") 
if (num.ln<=0)  stop ("num.ln must be a positive value") 

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
p <- ncol(data)-2
y <- data[, 1]
x <- data[, 2:(1+p)]
if (!is.matrix(x)) x <- as.matrix(x)
t <- data[, p+2]


if (is.null(w)) w <- c(quantile(data[,p+2],0.1), quantile(data[,p+2],0.9))



if ((b.equal.h==TRUE) & ((!is.null(num.b)) | (!is.null(num.h))) ) {num.b <- max(num.b, num.h); num.h <- num.b}
else if ((b.equal.h==TRUE) & (is.null(num.b)) & (is.null(num.h)) ) {num.b <- 50; num.h <- num.b}
else if ((b.equal.h==FALSE) & ((!is.null(num.b)) & (is.null(num.h))) ) num.h <- num.b
else if ((b.equal.h==FALSE) & ((is.null(num.b)) & (!is.null(num.h))) ) num.b <- num.h
else if ((b.equal.h==FALSE) & (is.null(num.b)) & (is.null(num.h)) ) {num.b <- 50; num.h <- num.b}


if ( is.null(b.seq) && is.null(h.seq) ) {
    a <- as.matrix(abs(outer(data[,p+2], data[,p+2],"-")))
    for (i in 1:n) {a[i,i] <- -1000}
    a <- as.vector(a[a!=-1000])

    b.min <- quantile(a,0.05)
    b.max <- (max(t)-min(t))*0.25

    b.seq <- seq(b.min,b.max,length.out=num.b)
  
    h.seq <- seq(b.min,b.max,length.out=num.h)
}

else if ( !is.null(b.seq) && !is.null(h.seq) && !b.equal.h) {
  num.b <- length(b.seq)
  num.h <- length(h.seq)
} 

else if ( !is.null(b.seq) && !is.null(h.seq) && b.equal.h) {
  c <- b.seq==h.seq
  if (any(c==FALSE)) stop("The input arguments b.seq and h.seq are not equal")
  if ( length(b.seq) != length(h.seq) ) stop("The input arguments b.seq and h.seq have different lengths")
  num.b <- length(b.seq)
  num.h <- length(h.seq)
} 

else if ( is.null(b.seq) && !is.null(h.seq) && b.equal.h ) {
  b.seq <- h.seq
  num.b <- length(h.seq)
  num.h <- length(h.seq)
} 

else if ( !is.null(b.seq) && is.null(h.seq) && b.equal.h ) {
  h.seq <- b.seq
  num.b <- length(b.seq)
  num.h <- length(b.seq)
} 

else if ( is.null(b.seq) && !is.null(h.seq) && !b.equal.h ) {
  num.h <- length(h.seq)
  
  a <- as.matrix(abs(outer(data[,p+2], data[,p+2],"-")))
  for (i in 1:n) {a[i,i] <- -1000}
  a <- as.vector(a[a!=-1000])
  
  b.min <- quantile(a,0.05)
  b.max <- (max(t)-min(t))*0.25
  
  b.seq <- seq(b.min,b.max,length.out=num.b)
} 

else if ( !is.null(b.seq) && is.null(h.seq) && !b.equal.h ) {
  num.b <- length(b.seq)
  
  a <- as.matrix(abs(outer(data[,p+2], data[,p+2],"-")))
  for (i in 1:n) {a[i,i] <- -1000}
  a <- as.vector(a[a!=-1000])
  
  b.min <- quantile(a,0.05)
  b.max <- (max(t)-min(t))*0.25
  
  h.seq <- seq(b.min,b.max,length.out=num.h)
} 



if (b.equal.h==TRUE) num.hh <- 1
else num.hh <- num.h



CV <- array(0,c(num.b, num.hh, num.ln))
CV.opt <- 1:num.ln
b.h.CV <- data.frame(matrix(0,3,num.ln),row.names=c("ln","b","h"))
  
ele.seq <- seq(from = ln.0, to =ln.0+(num.ln-1)*step.ln , by = step.ln)



BETAmat <- plrm.beta(data=data, b.seq=b.seq, estimator=estimator, kernel=kernel)$BETA


if (!is.matrix(BETAmat))  BETAmat <- t(as.matrix(BETAmat))
y <- as.vector(y)
YXB.mat <- y -x%*%BETAmat



if (estimator=="NW")
    
   for (k in 1:num.ln) {
      
      for (i in 1:n) {
        
         if ((w[1]<=t[i]) & (t[i]<=w[2])) {
          
            diff <- t-t[i]                 
            if (b.equal.h==TRUE) {Zmat <- matrix(rep(diff, num.b), nrow = num.b, byrow = T)}
            else {Zmat <- matrix(rep(diff, num.hh), nrow = num.hh, byrow = T)}
          
            Umat <- Zmat/h.seq
            Kmat <- kernel.function(Umat)
            Kmat[(Kmat<0) | abs(col(Kmat)-i)<=ele.seq[k]] <- 0
          
            S0 <- apply(Kmat, 1, sum)
          
            Kmat <- Kmat/S0

            if (b.equal.h==TRUE) {Kmat <- array(Kmat, c(num.b, n, num.hh)); Kmat <- aperm(Kmat, c(3,2,1))}
            else Kmat <- array(Kmat, c(num.hh, n, num.b))
                 
            Ymat <- t(YXB.mat)
            Ymat <- array(Ymat, c(num.b, n, num.hh));Ymat <- aperm(Ymat, c(3,2,1))
                    
            mhat <- t(apply(Ymat * Kmat, c(1,3), sum))
          
            regr <- matrix(rep(x[i,]%*%BETAmat, num.hh), nrow=num.b, byrow=F) + mhat
          
            CV[, , k] <- CV[, , k] + (y[i]-regr)^2
          
        } # if
     } # for
  } # for
  
  
else if (estimator=="LLP")
    
  for (k in 1:num.ln) {
      
     for (i in 1:n) {
        
        if ((w[1]<=t[i]) & (t[i]<=w[2])) {
          
           diff <- t-t[i]               
           if (b.equal.h==TRUE) {Zmat <- matrix(rep(diff, num.b), nrow = num.b, byrow = T)}
           else {Zmat <- matrix(rep(diff, num.hh), nrow = num.hh, byrow = T)}
           
           Umat <- Zmat/h.seq
           Kmat <- kernel.function(Umat)
           Kmat[(Kmat<0) | abs(col(Kmat)-i)<=ele.seq[k]] <- 0
          
           S0 <- apply(Kmat, 1, sum)
           S1 <- apply(Kmat*Zmat, 1, sum)
           S2 <- apply(Kmat*Zmat^2, 1, sum)
          
           Kmat <- Kmat * (S2 - Zmat*S1)/(S0*S2 - S1^2)

           if (b.equal.h==TRUE) {Kmat <- array(Kmat, c(num.b, n, num.hh));Kmat <- aperm(Kmat, c(3,2,1))}
           else {Kmat <- array(Kmat, c(num.hh, n, num.b))}
          
           Ymat <- t(YXB.mat)
           Ymat <- array(Ymat, c(num.b, n, num.hh));Ymat <- aperm(Ymat, c(3,2,1))
                   
           mhat <- t(apply(Ymat * Kmat, c(1,3), sum))  
          
           regr <- matrix(rep(x[i,]%*%BETAmat, num.hh), nrow=num.b, byrow=F) + mhat
          
           CV[, , k] <- CV[, , k] +(y[i]-regr)^2
          
        } # if
     } # i
  } # k
  
  
CV<-CV/n

for (k in 1:num.ln) {
  
    index.CV <- order(CV[, ,k])[1]
    
    if (b.equal.h==FALSE) {
    index.h <- 1+trunc((index.CV-1)/num.b);
    index.b <- index.CV - (index.h-1)*num.b}
        
    if (b.equal.h==TRUE) {b.h.CV[,k] <- c(ele.seq[k], b.seq[index.CV], h.seq[index.CV]); CV.opt[k] <- CV[index.CV, , k]}
    else {b.h.CV[,k] <- c(ele.seq[k], b.seq[index.b], h.seq[index.h]); CV.opt[k] <- CV[index.b, index.h, k]}
    
}    
    

list(bh.opt=b.h.CV, CV.opt=CV.opt, CV=CV, b.seq=b.seq, h.seq=h.seq, w=w)

}

