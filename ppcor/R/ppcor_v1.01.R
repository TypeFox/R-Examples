# partial correlation
pcor <- function(x, method = c("pearson", "kendall", "spearman"))
{
	# correlation method
	method <- match.arg(method)

	# check the data
	if (is.data.frame(x)) 
      	x <- as.matrix(x)
    	if (!is.matrix(x)) 
        	stop("supply a matrix-like 'x'")
    	if (!(is.numeric(x) || is.logical(x))) 
        	stop("'x' must be numeric")
    	stopifnot(is.atomic(x))

	# sample number
	n <- dim(x)[1]
	
	# given variables' number
	gp <- dim(x)[2]-2

	# covariance matrix
	cvx <- cov(x,method=method)

	# inverse covariance matrix
	if(det(cvx) < .Machine$double.eps){
    warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
		icvx <- ginv(cvx)
	}else
		icvx <- solve(cvx)

	# partial correlation
	pcor <- -cov2cor(icvx)
	diag(pcor) <- 1
	
	# p-value
	if(method == "kendall"){
		statistic <- pcor/sqrt(2*(2*(n-gp)+5)/(9*(n-gp)*(n-1-gp)))
		p.value <- 2*pnorm(-abs(statistic))

	}else{
		statistic <- pcor*sqrt((n-2-gp)/(1-pcor^2))
 		p.value <- 2*pt(-abs(statistic),(n-2-gp))
  		#p.value <- 2*pnorm(-abs(statistic))
	}

	diag(statistic) <- 0
	diag(p.value) <- 0

	list(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gp=gp,method=method)
}

# semi-partial (part) correlation
spcor <- function(x, method = c("pearson", "kendall", "spearman"))
{
	# correlation method
	method <- match.arg(method)

	# check the data
	if (is.data.frame(x)) 
      	x <- as.matrix(x)
    	if (!is.matrix(x)) 
        	stop("supply a matrix-like 'x'")
    	if (!(is.numeric(x) || is.logical(x))) 
        	stop("'x' must be numeric")
    	stopifnot(is.atomic(x))

	# sample number
	n <- dim(x)[1]
	
	# given variables' number
	gp <- dim(x)[2]-2

	# covariance matrix
	cvx <- cov(x,method=method)

	# inverse covariance matrix
	if(det(cvx) < .Machine$double.eps){
	  warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
	  icvx <- ginv(cvx)
	}else
	  icvx <- solve(cvx)

	# semi-partial correaltion
	spcor <- -cov2cor(icvx)/sqrt(diag(cvx))/sqrt(abs(diag(icvx)-t(t(icvx^2)/diag(icvx))))
	diag(spcor) <- 1

	# p-value
	if(method == "kendall"){
		statistic <- spcor/sqrt(2*(2*(n-gp)+5)/(9*(n-gp)*(n-1-gp)))
		p.value <- 2*pnorm(-abs(statistic))

	}else{
		statistic <- spcor*sqrt((n-2-gp)/(1-spcor^2))
 		p.value <- 2*pt(-abs(statistic),(n-2-gp))
  		#p.value <- 2*pnorm(-abs(statistic))
	}

	diag(statistic) <- 0
	diag(p.value) <- 0

	list(estimate=spcor,p.value=p.value,statistic=statistic,n=n,gp=gp,method=method)
}

# pairwise partial correlation
pcor.test <- function(x,y,z,method=c("pearson", "kendall", "spearman"))
{
	# The partial correlation coefficient between x and y given z
	#
	# pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
	#
	# x and y should be vectors
	#
	# z can be either a vector or a matrix

	# correlation method
	method <- match.arg(method)

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	# merge into a matrix
	xyz <- data.frame(x,y,z)

	# partial correlation
	pcor = pcor(xyz,method=method)

	data.frame(estimate=pcor$est[1,2],p.value=pcor$p.value[1,2],statistic=pcor$statistic[1,2],n=pcor$n,gp=pcor$gp,Method=method)
}	

# pairwise semi-partial (part) correlation
spcor.test <- function(x,y,z,method=c("pearson", "kendall", "spearman"))
{
	# The semi-partial (part) correlation coefficient between x and y given z
	#
	# spcor.test is free and comes with ABSOLUTELY NO WARRANTY.
	#
	# x and y should be vectors
	#
	# z can be either a vector or a matrix

	# correlation method
	method <- match.arg(method)

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	# merge into a matrix
	xyz <- data.frame(x,y,z)

	# semi-partial (part) correlation
	spcor = spcor(xyz,method=method)

	data.frame(estimate=spcor$est[1,2],p.value=spcor$p.value[1,2],statistic=spcor$statistic[1,2],n=spcor$n,gp=spcor$gp,Method=method)
}	
