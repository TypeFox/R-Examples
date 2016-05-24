
plrm.ancova <- function(data=data, t=t, b.seq=NULL, h.seq=NULL, w=NULL, estimator="NW", 
                        kernel="quadratic", time.series=FALSE, Var.Cov.eps=NULL, 
                        Tau.eps=NULL, b0=NULL, h0=NULL, lag.max=50, p.max=3, q.max=3, 
                        ic="BIC", num.lb=10, alpha=0.05)
{
if (!is.array(data))  stop("data must be an array")
if (ncol(data) < 2)  stop("data must have at least 2 columns")

if (is.null(t)) stop("t must not be NULL")
if (sum(is.na(t))  != 0) stop("t must have numeric values")
if (any(t<0)) stop ("all elements in t must be positive")
if (length(t)!= nrow(data)) stop ("length(t) and nrow(data) must be equal")

if ( (!is.null(b.seq)) && (sum(is.na(b.seq)) != 0) ) stop ("b.seq must be numeric")
if ( (!is.null(b.seq)) && (any(b.seq<=0)) ) stop ("b.seq must contain one ore more positive values")

if ( (!is.null(h.seq)) && (sum(is.na(h.seq)) != 0) ) stop ("h.seq must be numeric")
if ( (!is.null(h.seq)) && (any(h.seq<=0)) ) stop ("h.seq must contain one ore more positive values")

if ( (!is.null(w)) && (sum(is.na(w) )  != 0) ) stop ("w must be numeric")
if ( (!is.null(w)) && (!is.vector(w)) ) stop ("w must be a vector")
if ( (!is.null(w)) && (length(w)!=2) ) stop("w must be a vector of length 2")
if ( (!is.null(w)) && (any(w<0)) ) stop ("w must contain two positive values")
if ( (!is.null(w)) && (w[1]>w[2]) ) stop("w[2] must be greater than w[1]")
  
if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")

if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")

if (!is.logical(time.series)) stop("time.series must be logical")

if ( (!is.null(Var.Cov.eps)) && (sum(is.na(Var.Cov.eps)) != 0) ) stop("Var.Cov.eps must have numeric values")
if ( (!is.null(Var.Cov.eps)) && (!is.array(Var.Cov.eps)) )  stop("Var.Cov.eps must be an array")
if ( (!is.null(Var.Cov.eps)) && ( (ncol(Var.Cov.eps) != nrow(data)) | (nrow(Var.Cov.eps) != nrow(data)) | (dim(Var.Cov.eps)[3] != dim(data)[3]) ) ) stop("Var.Cov.eps must have dimension n x n x L")
if (!is.null(Var.Cov.eps)) {for (k in 1:(dim(Var.Cov.eps)[3]) ) if (any(t(Var.Cov.eps[,,k]) != Var.Cov.eps[,,k])) stop("Var.Cov.eps must be symmetric") }

if ( (!is.null(Tau.eps)) && (sum(is.na(Tau.eps))  != 0) ) stop("Tau.eps must have numeric values")
if ( (!is.null(Tau.eps)) && (!is.vector(Tau.eps)) )  stop("Tau.eps must be a vector")
if ( (!is.null(Tau.eps)) && (length(Tau.eps) != (ncol(data)-1)) ) stop ("Tau.eps must have length equals to ncol(data)-1")

if ( (!is.null(b0)) && (!is.numeric(b0)) )  stop ("b0 must be numeric")
if ( (!is.null(b0)) && (length(b0) !=1) ) stop ("b0 must be an only value") 
if ( (!is.null(b0)) && (b0<=0) ) stop ("b0 must be a positive value") 

if ( (!is.null(h0)) && (!is.numeric(h0)) )  stop ("h0 must be numeric")
if ( (!is.null(h0)) && (length(h0) !=1) ) stop ("h0 must be an only value") 
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
p <- ncol(data)-1
L <- dim(data)[3]

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
  
}
  
  

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


if ((is.null(Var.Cov.eps)) | (is.null(Tau.eps))) {

# We build auxiliar residuals in order to estimate Var.Cov.eps and/or Tau.eps
	eps.0 <- matrix(0,n,L)

      for (k in 1:L) {
   
      	y.x.t <- cbind(data[,,k],t)
      
     		beta.est.0 <- plrm.beta(data=y.x.t, b.seq=b0, estimator=estimator, kernel=kernel)$BETA

      	y.0 <- data[,1,k]-data[,-1,k]%*%beta.est.0
      
        m.0 <- np.est(data=cbind(y.0,t),newt=t, h.seq=h0, estimator=estimator, kernel=kernel)
         
		    eps.0[,k] <- y.0 - m.0
      
   	} # for
   
}
  
if (is.null(Var.Cov.eps)) {
  
    v.c.eps <- TRUE
   	
		Var.Cov.eps <- array(0,c(n,n,L))
		
		if (!time.series) {    

   			for (k in 1:L) {

   				var.eps <- var(eps.0[,k]) 
   				var.eps <- as.numeric(var.eps)	
  				Var.Cov.eps[,,k] <- diag(var.eps,n,n)
				
				} # for
		} # if
	
    else {    

   		for (k in 1:L) {
					V.eps <- matrix(NA, n, n)
					Var.Cov.mat <- var.cov.matrix(x=eps.0[,k], n=n, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)
					
					V.eps <- Var.Cov.mat[[1]]
					v.pv.Box.test <- Var.Cov.mat[[2]]
					v.pv.t.test <- Var.Cov.mat[[3]]
          v.ar.ma <- Var.Cov.mat[[4]]
          
					Var.Cov.eps[,,k] <- V.eps
			} # for
		} # else
} # if
else v.c.eps <- FALSE

if (is.null(Tau.eps)) {
  
    t.eps <- TRUE
  		
		if (!time.series) Tau.eps <- apply(eps.0, 2, var)  
	
    	else {    
				Var.Cov.sum <- var.cov.sum(X=eps.0, lag.max=lag.max, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)
				       
				Tau.eps <- Var.Cov.sum[[1]]
				t.pv.Box.test <- Var.Cov.sum[[2]]
				t.pv.t.test <- Var.Cov.sum[[3]]
        t.ar.ma <- Var.Cov.sum[[4]]
				
			} # else
} # if
else t.eps <- FALSE

#####################################
# PARAMETRIC AND NON-PARAMETRIC TEST
#####################################
YY <- matrix(0,n,L)
WY <- matrix(0,n,1)
  
XX <- array(0,c(n,p,L))
WX <- matrix(0,n,L)
  
data_par1 <- array(0,c(n,p+1,L))
  
Y <- matrix(0,n,L)
beta.est <- matrix(0,p,1)
  
  
for (i in 1:num.b) {
    
   for (k in 1:L) {
  
      data1 <- cbind(data[,1,k],t)
      
      WY <- np.est(data=data1,newt=t,h.seq=b.seq[i],estimator=estimator,kernel=kernel)
      YY[,k] <- data[,1,k]-WY
      data_par1[,1,k] <- YY[,k]
      
      
      y.x.t <- cbind(data[,,k],t)
      
      beta.est <- plrm.beta(data=y.x.t, b.seq=b.seq[i], estimator=estimator, kernel=kernel)$BETA
      Y[,k] <- data[,1,k]-data[,-1,k]%*%beta.est

      
      for (j in 1:p) {
         data2 <- cbind(data[,j+1,k],t)
        
         WX <- np.est(data=data2,newt=t,h.seq=b.seq[i],estimator=estimator,kernel=kernel)
         XX[,j,k] <- data[,j+1,k]-WX
         data_par1[,j+1,k] <- XX[,j,k]
      } # for
   } # for
    

   
# We test BETA_1= ... =BETA_L
par.test <- par.ancova(data=data_par1, time.series=time.series, Var.Cov.eps=Var.Cov.eps, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb) 

parametric.test[i,]  <- c(b.seq[i], par.test$par.ancova$Q.beta, par.test$par.ancova$p.value)
    
    
# We test m_1=m_2=...=m_L
data_nopar2 <- matrix(c(Y,t),nrow=n)
    
np.test <- np.ancova(data=data_nopar2, h0=h0, h.seq=h.seq[i], w=w, estimator=estimator, kernel=kernel,
                     time.series=time.series, Tau.eps=Tau.eps, lag.max=lag.max, p.max=p.max, q.max=q.max, ic=ic, alpha=alpha, num.lb=num.lb)

nonparametric.test[i,]  <- c(b.seq[i], h.seq[i], np.test$np.ancova$Q.m, np.test$np.ancova$Q.m.normalised, np.test$np.ancova$p.value)
        
} # for


parametric.test <- as.data.frame(parametric.test)
nonparametric.test <- as.data.frame(nonparametric.test)

  
if ( (v.c.eps) && (t.eps) && (time.series) ) list(parametric.test=parametric.test, nonparametric.test=nonparametric.test, pv.Box.test=v.pv.Box.test, pv.t.test=v.pv.t.test, v.ar.ma=v.ar.ma)
else if ( !(v.c.eps) && (t.eps) && (time.series) ) list(parametric.test=parametric.test, nonparametric.test=nonparametric.test, pv.Box.test=t.pv.Box.test, pv.t.test=t.pv.t.test, ar.ma=t.ar.ma)
else if ( (v.c.eps) && !(t.eps) && (time.series) ) list(parametric.test=parametric.test, nonparametric.test=nonparametric.test, pv.Box.test=v.pv.Box.test, pv.t.test=v.pv.t.test, ar.ma=v.ar.ma)
else list(parametric.test=parametric.test, nonparametric.test=nonparametric.test)

}

