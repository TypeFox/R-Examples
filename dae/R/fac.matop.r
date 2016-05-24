"fac.ar1mat" <- function(factor, rho)
#function to form the correlation matrix for a (generalized) factor where the 
#  correlation between the levels has an ar1 pattern with correlation parameter 
#  rho;
#the order of the matrix (n) is the length of the factor and
#the power of the exponent of rho is derived from the difference in the levels
#of the factor, these levels being converted to numeric.
#The method 
#  a) forms an n x n matrix of all pairwise differences in the numeric values 
#     corresponding to the observed levels of the factor by taking the 
#     difference between the following two n x n matrices are equal: 1) each row 
#     contains the numeric values corresponding to the observed levels of the 
#     factor, and 2) each column contains the numeric values corresponding to 
#     the observed levels of the factor, 
#  b) replaces each element of the pairwise difference matrix with rho raised to 
#     the absolute value of the difference.
{ if (!is.factor(factor))
    stop("Must supply a single factor as the argument")
	repl <- table(factor)
  lev.fac <- levels(factor)
  if (0 %in% repl)
  { lev.fac <- names(repl)[!(repl == 0)]
  }
  fact.new <- factor(factor, labels=lev.fac)
  numfac <- as.numfac(fact.new)
  order <- length(fact.new)
  n <- order*order
  ar1 <- matrix(rep(rho, n), nrow=order, ncol=order)
  row.no <- matrix(rep(numfac, times=order), nrow=order, ncol=order)
  col.no <- matrix(rep(numfac, each=order), nrow=order, ncol=order)
  power <- abs(row.no - col.no)
  ar1 <- ar1^power
  ar1
}

"fac.meanop" <- function(factor)
{
#computes the projection matrix that produces the means corresponding to a (generalized) factor
#The method conputes the design matrix, X, for the factor and then X(diag(1/reps))t(X)
  if (!is.factor(factor))
    stop("Must supply a single factor as the argument")
  n <- length(factor)
	repl <- table(factor)
  lev.fac <- levels(factor)
  if (0 %in% repl)
  { lev.fac <- names(repl)[!(repl == 0)]
  }
  fact.new <- factor(factor, labels=lev.fac)
	l.fac <- length(lev.fac)
	X.fac <- matrix(rep(fact.new, times=l.fac), ncol=l.fac, dimnames=list(1:n, lev.fac))
	X.fac <- sapply(1:l.fac, function(i, X.fac) as.numeric(X.fac[,i] == lev.fac[i]), X.fac=X.fac)
	repl.new <- table(fact.new)
	div <- diag(1/repl.new, nrow=length(repl.new))
	M.fac <- X.fac %*% div %*% t(X.fac)
	M.fac <- projector(M.fac)
	if (degfree(M.fac) != l.fac)
	   stop("Degrees of freedom of projector not equal to observed number of levels of ",deparse(substitute(factor)),"\n")
	M.fac
}

"fac.sumop" <- function(factor)
{
#computes the summation matrix corresponding to a (generalized) factor
#The method conputes the design matrix, X, for the factor and then X t(X)
#An alternative method is that in fac.vcmat
  if (!is.factor(factor))
    stop("Must supply a single factor as the argument")
  n <- length(factor)
	repl <- table(factor)
  lev.fac <- levels(factor)
  if (0 %in% repl)
  { lev.fac <- names(repl)[!(repl == 0)]
  }
  fact.new <- factor(factor, labels=lev.fac)
	l.fac <- length(lev.fac)
	X.fac <- matrix(rep(fact.new, times=l.fac), ncol=l.fac, dimnames=list(1:n, lev.fac))
	X.fac <- sapply(1:l.fac, function(i, X.fac) as.numeric(X.fac[,i] == lev.fac[i]), X.fac=X.fac)
	S.fac <- X.fac %*% t(X.fac)
	S.fac
}

"fac.vcmat" <- function(factor, sigma2)
#function to form the variance matrix for the variance component of a 
#  (gener4alized) factor;
#that is, for a factor for which effects for different levels are independent 
#  and identically distributed with their variance given by the variance component;
#  elements of the matrix will equal either zero or sigma2 and displays 
#  compound symmetry.
#the order of the matrix (n) is the length of the factor and
#sigma2 is the value of the variance component.
#The method 
#  a) forms the n x n summation or relationship matrix whose elements are equal 
#     to zero except for those elements whose corresponding elements in the 
#     following two n x n matrices are equal: 1) each row contains the numeric 
#     values corresponding to the observed levels of the factor, and 2) each 
#     column contains the numeric values corresponding to the observed levels 
#     of the factor, 
#  b) multiplies the summation matrix by sigma2.
{ if (!is.factor(factor))
    stop("Must supply a single factor as the argument")
	repl <- table(factor)
  lev.fac <- levels(factor)
  if (0 %in% repl)
  { lev.fac <- names(repl)[!(repl == 0)]
  }
  fact.new <- factor(factor, labels=lev.fac)
  numfac <- as.numfac(fact.new)
  order <- length(fact.new)
  n <- order*order
  id <- matrix(rep(0, n), nrow=order, ncol=order)
  row.no <- matrix(rep(numfac, times=order), nrow=order, ncol=order)
  col.no <- matrix(rep(numfac, each=order), nrow=order, ncol=order)
  id[row.no == col.no] <- 1
  id <- sigma2*id
  id
}
