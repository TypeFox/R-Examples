
plrm.gof <- function(data=data, beta0=NULL, m0=NULL, b.seq=NULL, h.seq=NULL, 
                     w=NULL, estimator="NW", kernel="quadratic", time.series=FALSE, 
                     Var.Cov.eps=NULL, Tau.eps=NULL, b0=NULL, h0=NULL, lag.max=50, 
                     p.max=3, q.max=3, ic="BIC", num.lb=10, alpha=0.05) 
{

if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data)<3)  stop("data must have at least 3 columns: y, x, t")

if ( (!is.null(beta0)) && (sum(is.na(beta0))  != 0) ) stop("beta0 must have numeric values")
if ( (!is.null(beta0)) && (!is.vector(beta0)) ) stop("beta0 must be a vector")
if ( (!is.null(beta0)) && (length(beta0) != (ncol(data)-2)) ) stop("beta0 must have length equals to ncol(data)-2")

if ( (!is.null(m0)) && (!is.function(m0)) ) stop ("m0 must be a function")

if ( (!is.null(b.seq)) && (sum(is.na(b.seq))  != 0) ) stop ("b.seq must be numeric")
if ( (!is.null(b.seq)) && (any(b.seq<=0)) ) stop ("b.seq must contain one ore more positive values")

if ( (!is.null(h.seq)) && (sum(is.na(h.seq))  != 0) ) stop ("h.seq must be numeric")
if ( (!is.null(h.seq)) && (any(h.seq<=0)) ) stop ("h.seq must contain one ore more positive values")

if ( (!is.null(w)) && (sum(is.na(w) )  != 0) ) stop ("w must be numeric")
if ( (!is.null(w)) && (!is.vector(w)) ) stop ("w must be a vector")
if ( (!is.null(w)) && (length(w)!=2) ) stop("w must be a vector of length 2")
if ( (!is.null(w)) && (any(w<0)) ) stop ("w must contain two positive values")
if ( (!is.null(w)) && (w[1]>w[2]) ) stop("w[2] must be greater than w[1]")

if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")

if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")

if (!is.logical(time.series)) stop("time.series must be logical")

if ( (!is.null(Var.Cov.eps)) && (sum(is.na(Var.Cov.eps))  != 0) ) stop("Var.Cov.eps must have numeric values")
if ( (!is.null(Var.Cov.eps)) && (!is.matrix(Var.Cov.eps)) )  stop("Var.Cov.eps must be a matrix")
if ( (!is.null(Var.Cov.eps)) && ( (ncol(Var.Cov.eps) != nrow(data)) | (nrow(Var.Cov.eps) != nrow(data)) ) ) stop("Var.Cov.eps must have dimension n x n")
if ( (!is.null(Var.Cov.eps)) && (any(t(Var.Cov.eps) != Var.Cov.eps)  ) ) stop("Var.Cov.eps must be symmetric")

if ( (!is.null(Tau.eps)) && (length(Tau.eps) !=1) ) stop ("Tau.eps must be an only value")
if ( (!is.null(Tau.eps)) && (!is.numeric(Tau.eps)) )  stop ("Tau.eps must be numeric") 

if ( (!is.null(b0)) && (length(b0) !=1) ) stop ("b0 must be an only value") 
if ( (!is.null(b0)) && (!is.numeric(b0)) )  stop ("b0 must be numeric")
if ( (!is.null(b0)) && (b0<=0) ) stop ("b0 must be a positive value") 

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
p <- ncol(data)-2
t <- data[,p+2]

if (is.null(beta0)) beta0 <- rep(0,length.out=p)

if (is.null(m0)) {m0 <- function(u) {0} }

if ( (is.null(h.seq)) & (!is.null(b.seq)) ) {h.seq <- b.seq}
if ( (is.null(b.seq)) & (!is.null(h.seq)) ) {b.seq <- h.seq}
if (length(b.seq) != length(h.seq)) stop("length(b.seq) and length(h.seq) must be equal")

if ( (is.null(h0)) & (!is.null(b0)) ) {h0 <- b0}
if ( (is.null(b0)) & (!is.null(h0)) ) {b0 <- h0}


  
dif <- t[2]-t[1]
                       
for (i in 1:(n-1)) {
     
  dist <- t[i+1]-t[i]
  if (abs(dist - dif) < 2e-15) {}
  else stop("t values must be equidistant")    
      
} # for
    


x0 <- min(t)
x1 <- max(t)

y0 <- (1-0.5)/n
y1 <- (n-0.5)/n

slope <- (y1 - y0)/(x1-x0)
intercept <- y1-slope*x1
  


if (is.null(b.seq)) {b.seq <- (seq(0.05, 0.25, length.out=10))/slope; h.seq <- b.seq}
num.h <- length(h.seq)
num.b <- length(b.seq)

if (is.null(b0)) {b0 <- 0.25/slope; h0 <- b0}

if (is.null(w)) w <- (-intercept + c(0.1, 0.9))/slope



parametric.test <- matrix(1000,num.b,3)
dimnames(parametric.test) <- list(NULL, c("b", "Q", "p.v"))
    
nonparametric.test <- matrix(1000,num.h,5)
dimnames(nonparametric.test) <- list(NULL,c("b", "h", "Q", "Q.normalised", "p.v"))
  
  
#####################################
# PARAMETRIC AND NON-PARAMETRIC TEST
#####################################
YY <- matrix(0,n,1)
WY <- matrix(0,n,1)
  
XX <- matrix(0,n,p)
WX <- matrix(0,n,1)
  
data_par <- matrix(0,n,p+1)
data_nopar <- matrix(0,n,2)
  
Y <- matrix(0,n,1)
beta.est <- matrix(0,p,1)
  
  

if ((is.null(Var.Cov.eps)) | (is.null(Tau.eps))) {

# We build auxiliar residuals in order to estimate Var.Cov.eps and/or Tau.eps
         
    beta.est.0 <- plrm.beta(data=data, b.seq=b0, estimator=estimator, kernel=kernel)$BETA

    # y.0 <- data[,1]-data[,-c(1,p+2)]%*%beta.est.0
    y.0 <- data[,1]-as.matrix(data[,-c(1,p+2)])%*%beta.est.0
    
      
    m.0 <- np.est(data=cbind(y.0,t),newt=t, h.seq=h0, estimator=estimator, kernel=kernel)
         
		eps.0 <- y.0 - m.0
    
} # if



if (is.null(Var.Cov.eps)) {
    v.c.eps <- TRUE
		
	  if (!time.series) {    

   		  var.eps <- var(eps.0) 
   			var.eps <- as.numeric(var.eps)	
  			Var.Cov.eps <- diag(var.eps,n,n)
										
		}	# if
    else {    

				Var.Cov.eps <- matrix(NA, n, n)
				Var.Cov.mat <- var.cov.matrix(x=eps.0, n=n, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)
				
				Var.Cov.eps <- Var.Cov.mat[[1]]
				v.pv.Box.test <- Var.Cov.mat[[2]]
				v.pv.t.test <- Var.Cov.mat[[3]]
        v.ar.ma <- Var.Cov.mat[[4]]
        
		} #else
} # if
else v.c.eps <- FALSE

if (is.null(Tau.eps)) {
  		
    t.eps <- TRUE
		if (!time.series) Tau.eps <- var(eps.0)  
	
    else {  
      
				Var.Cov.sum <- var.cov.sum(X=eps.0, lag.max=lag.max, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)
				
				Tau.eps <- Var.Cov.sum[[1]]
				t.pv.Box.test <- Var.Cov.sum[[2]]
				t.pv.t.test <- Var.Cov.sum[[3]]
        t.ar.ma <- Var.Cov.sum[[4]]
    
		} # else
} # if
else t.eps <- FALSE

for (i in 1:num.b) {
    
    data1 <- cbind(data[,1],t)
      
    WY <- np.est(data=data1,newt=t,h.seq=b.seq[i],estimator=estimator,kernel=kernel)
    YY <- data[,1]-WY
      
      
    beta.est <- plrm.beta(data=data, b.seq=b.seq[i], estimator=estimator, kernel=kernel)$BETA
    # Y <- data[,1]-data[,2:(p+1)]%*%beta.est
    Y <- data[,1]-as.matrix(data[,2:(p+1)])%*%beta.est
    

      
    for (j in 1:p) {
      data2 <- cbind(data[,j+1],t)
        
      WX <- np.est(data=data2,newt=t,h.seq=b.seq[i],estimator=estimator,kernel=kernel)
      XX[,j] <- data[,j+1]-WX
    } # for
            
    
    
    # We test BETA = BETA0
    data_par <- matrix(c(YY,XX),ncol=p+1)   
            
    par.test <- par.gof(data=data_par, beta0=beta0, time.series=time.series, Var.Cov.eps=Var.Cov.eps, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)       
           
    parametric.test[i,]  <- c(b.seq[i], par.test$par.gof$Q.beta, par.test$par.gof$p.value)
    
    
    # We test m = m0
    data_nopar <- matrix(c(Y,t),ncol=2)
    np.test <- np.gof(data=data_nopar, m0=m0, h0=h0, h.seq=h.seq[i], w=w, estimator=estimator, kernel=kernel,
                       time.series=time.series, Tau.eps=Tau.eps, lag.max=lag.max, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)

    nonparametric.test[i,]  <- c(b.seq[i], h.seq[i], np.test$np.gof$Q.m, np.test$np.gof$Q.m.normalised, np.test$np.gof$p.value)
    
} # for

parametric.test <- as.data.frame(parametric.test)
nonparametric.test <- as.data.frame(nonparametric.test)


if ( (v.c.eps) && (t.eps) && (time.series) ) list(parametric.test=parametric.test, nonparametric.test=nonparametric.test, pv.Box.test=v.pv.Box.test, pv.t.test=v.pv.t.test, ar.ma=v.ar.ma)
else if ( !(v.c.eps) && (t.eps) && (time.series) )  list(parametric.test=parametric.test, nonparametric.test=nonparametric.test, pv.Box.test=t.pv.Box.test, pv.t.test=t.pv.t.test, ar.ma=t.ar.ma)
else if ( (v.c.eps) && !(t.eps) && (time.series) )  list(parametric.test=parametric.test, nonparametric.test=nonparametric.test, pv.Box.test=v.pv.Box.test, pv.t.test=v.pv.t.test, ar.ma=v.ar.ma)
else list(parametric.test=parametric.test, nonparametric.test=nonparametric.test)

}

