
plrm.beta <- function(data=data, b.seq=NULL, 
                      estimator="NW", kernel = "quadratic")
{

if (!is.matrix(data))  stop("data must be a matrix")
if (ncol(data)<3)  stop("data must have at least 3 columns: y, x, t")

if ( (!is.null(b.seq)) && (sum(is.na(b.seq))  != 0) ) stop ("b.seq must be numeric")
if ( (!is.null(b.seq)) && (any(b.seq<=0)) ) stop ("b.seq must contain one ore more positive values")

if ((estimator != "NW") & (estimator != "LLP"))  stop("estimator=NW or estimator=LLP is required")

if ((kernel != "quadratic") & (kernel != "Epanechnikov") & (kernel != "triweight") & (kernel != "gaussian") & (kernel != "uniform"))  stop("kernel must be one of the following: quadratic, Epanechnikov, triweight, gaussian or uniform")



n <- nrow(data)
p <- ncol(data)-2
y <- data[,1]
x <- data[,2:(p+1)]
if (!is.matrix(x)) x <- as.matrix(x)
t <- data[,p+2]


if (is.null(b.seq)) b.seq <- plrm.cv(data=data, estimator=estimator, kernel=kernel)$bh.opt[2,1]
num.b <- length(b.seq)

     
XX <- array(0,c(n,p,num.b))
G <- array(0,c(n,p,num.b))
BETA <- matrix(0,p,num.b)


for (j in 1:p) {
  
     G[,j,] <- np.est(data=data[,c(1+j,p+2)], newt=t, h.seq=b.seq, estimator=estimator, kernel=kernel)
     XX[,j,]<- x[,j]-G[,j,]
     
} # for


y <- as.vector(y)
YY <- y-np.est(data=data[,c(1,p+2)], newt=t, h.seq=b.seq, estimator=estimator, kernel=kernel)


for (k in 1:num.b) {
  
     if ((sum(YY[,k])==Inf) | (sum(YY[,k])=="NaN")) BETA[,k]<-0/0
            
     else {
          yy.xx <- cbind(YY[,k], XX[,,k])
          BETA[,k] <- par.est(data=yy.xx)
      } # else
     
} # for


list(BETA=BETA, G=G)

}

