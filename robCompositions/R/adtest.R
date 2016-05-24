#' Anderson-Darling Normality Tests
#' 
#' This function provides three kinds of Anderson-Darling Normality Tests
#' (Anderson and Darling, 1952).
#' 
#' Three version of the test are implemented (univariate, angle and radius
#' test) and it depends on the data which test is chosen.
#' 
#' If the data is univariate the univariate Anderson-Darling test for normality
#' is applied.
#' 
#' If the data is bivariate the angle Anderson-Darling test for normality is
#' performed out.
#' 
#' If the data is multivariate the radius Anderson-Darling test for normality
#' is used.
#' 
#' If \sQuote{locscatt} is equal to \dQuote{robust} then within the procedure,
#' robust estimates of mean and covariance are provided using \sQuote{covMcd()}
#' from package robustbase.
#' 
#' To provide estimates for the corresponding p-values, i.e. to compute the
#' probability of obtaining a result at least as extreme as the one that was
#' actually observed under the null hypothesis, we use Monte Carlo techniques
#' where we check how often the statistic of the underlying data is more
#' extreme than statistics obtained from simulated normal distributed data with
#' the same (column-wise-) mean(s) and (co)variance.
#' 
#' @param x either a numeric vector, or a data.frame, or a matrix
#' @param R Number of Monte Carlo simulations to obtain p-values
#' @param locscatt standard for classical estimates of mean and (co)variance.
#' robust for robust estimates using \sQuote{covMcd()} from package robustbase
#' @return \item{statistic }{The result of the corresponding test statistic}
#' \item{method }{The chosen method (univariate, angle or radius)}
#' \item{p.value }{p-value}
#' @note These functions are use by \code{\link{adtestWrapper}}.
#' @author Karel Hron, Matthias Templ
#' @seealso \code{\link{adtestWrapper}}
#' @references Anderson, T.W. and Darling, D.A. (1952) Asymptotic theory of
#' certain goodness-of-fit criteria based on stochastic processes. \emph{Annals
#' of Mathematical Statistics}, \bold{23} 193-212.
#' @keywords htest
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom robustbase covMcd
#' @examples
#' 
#' adtest(rnorm(100))
#' data(machineOperators)
#' x <- machineOperators
#' adtest(isomLR(x[,1:2]))
#' adtest(isomLR(x[,1:3]))
#' adtest(isomLR(x))
#' adtest(isomLR(x[,1:2]), locscatt="robust")
#' 
adtest=function(x, R=1000, locscatt="standard") {   
	#DNAME <- deparse(substitute(x))	
	if( R < 1 ) stop("choose a higher value for R")
	if( R < 50 ) warnings("maybe, the estimation of the p-value(s) is not accurate; choose a higher value for R")
	cv <- function(x, type) {
		classical <- function(x){
			if( (length(dim(x)) < 1) | is.vector(x) ){
	 			list(mean=mean(as.numeric(x), na.rm=TRUE), varmat=var(x, na.rm=TRUE))
			} else {
				list(mean=colMeans(x, na.rm=TRUE), varmat=cov(x))  
			}
		}

		robust <- function(x){
			if( (length(dim(x)) < 1) | is.vector(x) ){
				list(mean=median(x), varmat=mad(x)^2)
			} else {
				v <- robustbase::covMcd(x)
				list(mean=v$center, varmat=v$cov)
			}
		}
		switch(type,
				standard = classical(x),
				robust = robust(x))
	}
	
	centre <- function(x, type) {
		switch(type,
				mean = mean(as.numeric(x)),
				median = median(as.numeric(x)),
				trimmed = mean(as.numeric(x), trim = .1))
	}
    if(locscatt=="standard") location <- "mean" else location <- "median"    

	### 1-dim:
	if((length(dim(x)) < 1) | is.vector(x) ){
	  x <- as.vector(x)
	  n <- length(x)
      if (n < 8 ) 
        stop("sample size must be greater than 7")
      stat=function(x, N=n, location="mean"){
        x <- sort(x[complete.cases(x)])
		estCv <- cv(x, locscatt)
        p <- pnorm((x - estCv$mean)/sqrt(estCv$varmat))
        h <- (2 * seq(1:N) - 1) * (log(p) + log(1 - rev(p)))
        A <- (25/N^2-4/N-1)*(centre(h, location)+N)
      }
      A=stat(x, location=location)
	  estCv <- cv(x, locscatt)
	  mv <- estCv$mean
	  varmat <- estCv$varmat
      n=length(c(x))
	  p <- sapply(X=1:R, FUN=function(X,...){ stat(rnorm(n, mv, sqrt(varmat)), location=location) })
      pvalue=mean(p>A)
      RVAL <- list(statistic = c(A = A), method = "A-D univariate normality test", 
        p.value=pvalue)
    } else if (dim(x)[2] == 2){
	### 2-dim:
	  n <- nrow(x)
	  stat=function(x, N=n, location="mean"){
		estCv <- cv(x, locscatt)
		varmat <- estCv$varmat
		mu <- estCv$mean
		u <- (1/sqrt(det(varmat)))*((x[,1]-mu[1])*sqrt(varmat[2,2])-(x[,2]-mu[2])*(varmat[1,2]/sqrt(varmat[2,2])))
		v <- (x[,2]-mu[2])/sqrt(varmat[2,2])
		teta <- atan(v/u)+(1-sign(u))*pi/2+(1+sign(u))*(1-sign(v))*pi/2
		z=teta/(2*pi)
		p=sort(z)
		h <- (2 * seq(1:N) - 1) * (log(p) + log(1 - rev(p)))
		A <- -N-centre(h, location)
	  }
	A=stat(x, location=location)  
	estCv <- cv(x, locscatt)
	varmat <- estCv$varmat
	mv <- estCv$mean	
	p <- sapply(X=1:R, FUN=function(X,...){ stat(mvrnorm(n, mv, sqrt(varmat)), location=location) })
	pvalue=mean(p>A)
	RVAL <- list(statistic = c(A = A), method = "A-D bivariate normality test", 
			 p.value=pvalue)
  } else { 
	### >= 3:
		n <- nrow(x)		
		stat=function(x, N=n, location="mean"){ 
			estCv <- cv(x, locscatt)
			#varmat <- var(x)
			#mu <- apply(x,2,mean)
			u <- mahalanobis(x, center=estCv$mean, cov=estCv$varmat)
			z <- pchisq(u,ncol(x))
			p=sort(z) 
			h <- (2 * seq(1:N) - 1) * (log(p) + log(1 - rev(p)))
			A <- -N-centre(h, location)
#par(mfrow=c(1,2)); plot(p) ; plot(h)
		} 
		A=stat(x, location=location)    
		#print(paste("A =", A))
		estCv <- cv(x, locscatt)		
		# estCv$mean=colMeans(x)  #weg
		 #estCv$varmat=var(x) #weg
		p=numeric(R)
		n=nrow(x)
		p <- sapply(X=1:R, FUN=function(X,...){ stat(mvrnorm(n, estCv$mean, estCv$varmat), location=location) })
		pvalue=mean(p>A)
		RVAL <- list(statistic = A, method = "A-D radius test", 
				p.value=pvalue)
  }		
		
	class(RVAL) <- "htest"
	return(RVAL)
	
}


######################################################################
#Cramer-vom Mises
#
#cvmtest=function (x) 
#{
#    DNAME <- deparse(substitute(x))
#    x <- sort(x[complete.cases(x)])
#    n <- length(x)
#    if (n < 8) 
#        stop("sample size must be greater than 7")
#    p <- pnorm((x - mean(x))/sd(x))
#    W <- (1/(12 * n) + sum((p - (2 * seq(1:n) - 1)/(2 * n))^2))*((2*n+1)/(2*n))
#    RVAL <- list(statistic = c(W = W), method = "Cramer-von Mises normality test", 
#        data.name = DNAME)
#    class(RVAL) <- "htest"
#    return(RVAL)
#}
#
#######################################################################
##Watson test
#
#wattest=function (x) 
#{
#    DNAME <- deparse(substitute(x))
#    x <- sort(x[complete.cases(x)])
#    n <- length(x)
#    if (n < 8) 
#        stop("sample size must be greater than 7")
#    p <- pnorm((x - mean(x))/sd(x))
#    W <- (1/(12 * n) + sum((p - (2 * seq(1:n) - 1)/(2 * n))^2))*((2*n+1)/(2*n))
#    WW <- W-((2*n+1)/2)*(mean(p)-1/2)^2
#    RVAL <- list(statistic = c(WW = WW), method = "Watson normality test", 
#        data.name = DNAME)
#    class(RVAL) <- "htest"
#    return(RVAL)
#}

