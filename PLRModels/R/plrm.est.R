
plrm.est <- function(data=data, b=NULL, h=NULL, newt=NULL, 
                     estimator="NW", kernel = "quadratic")
{

if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data)<3)  stop("data must have at least 3 columns: y, x, t")

if ( (!is.null(b)) && (length(b) !=1) ) stop ("b must be an only value") 
if ( (!is.null(b)) && (!is.numeric(b)) )  stop ("b must be numeric")
if ( (!is.null(b)) && (b<=0) ) stop ("b must be a positive value") 

if ( (!is.null(h)) && (length(h) !=1) ) stop ("h must be an only value") 
if ( (!is.null(h)) && (!is.numeric(h)) )  stop ("h must be numeric")
if ( (!is.null(h)) && (h<=0) ) stop ("h must be a positive value") 
  
if ( (!is.null(newt)) && (sum(is.na(newt))  != 0) ) stop ("newt must be numeric")
if ( (!is.null(newt)) && (any(newt<=0)) ) stop ("newt must contain one ore more positive values")

if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")

if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")


  
if ( (is.null(b)) & (is.null(h)) ) {
  b <- plrm.cv(data=data, estimator=estimator, kernel=kernel)$bh.opt[2,1]    
  h <- b}

else if ( (!is.null(b)) & (is.null(h)) ) h <- b
else if ( (is.null(b)) & (!is.null(h)) ) b <- h



beta <- plrm.beta(data=data, b.seq=b, estimator=estimator, kernel=kernel)$BETA
  
n <- nrow(data)
p <- ncol(data)-2
y <- data[, 1]
x <- data[, 2:(p+1)]
if (!is.matrix(x)) x <- as.matrix(x)
t <- data[, p+2]
  


if (!is.null(newt)) {m.newt <- rep(0, length(newt))}
m.t <- rep(0,length=n)
  

y <- as.vector(y)
data2 <- y-x%*%beta
data3 <- cbind(data2,t)
if (!is.null(newt)) {m.newt <- np.est(data=data3, h.seq=h, newt=newt, estimator=estimator, kernel=kernel)}
m.t <- np.est(data=data3, h.seq=h, newt=t, estimator=estimator, kernel=kernel) 
    
  

res <- y-x%*%beta-m.t
      
fitted.values <- x%*%beta+m.t


if (!is.null(newt)) {list(beta=beta, m.t=m.t, m.newt=m.newt, residuals=res, fitted.values=fitted.values, b=b, h=h)}
else {list(beta=beta, m.t=m.t, residuals=res, fitted.values=fitted.values, b=b, h=h)}
  
}

