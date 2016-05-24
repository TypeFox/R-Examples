#' @title Measure an impact of the covariates on the response using Lasso
#' This function evaluates the Lasso coefficients regressing \code{y} onto the design matrix \code{x} over subsamples in \code{subsamples}.
#' @details To solve the Lasso problem, we implement the coordinate descent algorithm as in Breheny Jian (2011). 
#' @param x Matrix with \code{n} observations of \code{p} covariates in each row.
#' @param y Response vector with \code{n} observations.
#' @param subsamples Matrix with \code{m} indices of \code{N} subsamples in each column. 
#' @param nonzero Number of non-zero coefficients estimated for each subsample.
#' @param family Determines the likelihood optimised in the estimation procedure.
#' @param alpha Scalar between 0 and 1 determining l2 penalty (see details).
#' @param maxit Maximum number of itarations when computing the lasso coefficients.
#' @param tol Scalar determining convergence of the lasso algorithm used.
#' @param lambda.ratio Scalar being a fraction of 1. Used in the lasso algorithm
#' @param nlam Number of penalty parameters used in the lasso algorithm.
#' @param ... Not in use.
#' @useDynLib rbvs lasso_coef_gaussian_r lasso_coef_binomial_r
#' @export 
#' @author Rafal Baranowski, Patrick Breheny
#' @references Tibshirani, Robert. "Regression shrinkage and selection via the lasso." Journal of the Royal Statistical Society. Series B (Methodological) (1996): 267-288.
#' @references Breheny, Patrick, and Jian Huang. "Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection." The Annals of Applied Statistics 5.1 (2011): 232.


lasso.coef <- function(x,y,subsamples,
                              nonzero=NULL,
                              family=c("gaussian", "binomial"),
                              alpha=1.0,
                              maxit=500,
                              tol=1e-2,
                              lambda.ratio=1e-6,
                              nlam=25,
                              ...){

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  y <- as.numeric(y)
  storage.mode(y) <- "double"
  
  subsamples <- as.matrix(subsamples)
  storage.mode(subsamples) <- "integer"
  
  
  #check if the sizes of x and y match
  
  n <- length(y)
  if(nrow(x) != n) stop("The number of rows in 'x' must be equal to the length of 'y'.")
  
  #check for NA's
  if(any(is.na(x))) stop("Matix 'x' cannot contain NA's.")
  if(any(is.na(y))) stop("Vector 'y' cannot contain NA's.")
  if(any(is.na(subsamples))) stop("Matrix 'subsamples' cannot contain NA's.")

  alpha <- as.double(alpha)[1]
  
  if(alpha < 0 || alpha > 1) stop("Parameter 'alpha' has to be between 0 and 1.")
  
  maxit <- as.integer(maxit)[1]

  if(maxit < 50){
    warning("Paramater 'maxit' must be greater than or equal to 50. Setting 'maxit' to 50.")
    maxit <- 50
  } 
  
  tol <- as.double(tol)[1]
  
  if(tol <= 0) stop("'tol' has to be positive.") 
  
  
  lambda.ratio <- as.double(lambda.ratio)[1]
  if(lambda.ratio < 0 || lambda.ratio > 1) stop("'lambda.ratio' has to be between 0 and 1.")

  nlam <- as.integer(nlam)[1]
  
  if(nlam < 2) {
    warning("Paramater 'nlam' must be greater than or equal to 50. Setting 'nlam' to 2.")
    maxit <- 2
  }

  
  #check if subsamples are actually subsamples
  min.ind <- min(subsamples)
  max.ind <- max(subsamples)
  if(min.ind < 1 || max.ind > n) stop("Elements of 'subsamples' must be between 1 and length(y).")
  
  

  
  m <- nrow(subsamples);
  B <- ncol(subsamples);
  p <- ncol(x)
  n <- nrow(x)
  
  if(is.null(nonzero)) nonzero <- m-1;
  

 
  
  family <- match.arg(family)
  
  
  if( family=="gaussian"){
    return(.Call("lasso_coef_gaussian_r", 
                 as.matrix(subsamples),
                 as.matrix(x),
                 as.numeric(y),
                 as.integer(nonzero),
                 as.double(alpha),
                 as.double(tol),
                 as.integer(maxit),
                 as.double(lambda.ratio),
                 as.integer(nlam)))   
  }else{
    return(.Call("lasso_coef_binomial_r", 
                 as.matrix(subsamples),
                 as.matrix(x),
                 as.numeric(y),
                 as.integer(nonzero),
                 as.double(alpha),
                 as.double(tol),
                 as.integer(maxit),
                 as.double(lambda.ratio),
                 as.integer(nlam)))
  }
  

  
  
}


#' @title Measure an impact of the covariates on the response using MC+.
#' This function evaluates the MC+ coefficients regressing \code{y} onto the design matrix \code{x} over subsamples in \code{subsamples}.
#' @details To solve the MC+ problem, we implement the coordinate descent algorithm as in Breheny Jian (2011). 
#' @param x Matrix with \code{n} observations of \code{p} covariates in each row.
#' @param y Response vector with \code{n} observations.
#' @param subsamples Matrix with \code{m} indices of \code{N} subsamples in each column. 
#' @param nonzero Number of non-zero coefficients estimated for each subsample.
#' @param family Determines the likelihood optimised in the estimation procedure.
#' @param alpha Scalar between 0 and 1 determining l2 penalty (see details).
#' @param gamma Scalar greater than 1. The concacivity parameter (see details).
#' @param maxit Maximum number of itarations when computing the MC+ coefficients.
#' @param tol Scalar determining convergence of the MC+ algorithm used.
#' @param lambda.ratio Scalar being a fraction of 1. Used in the MC+ algorithm
#' @param nlam Number of penalty parameters used in the MC+ algorithm.
#' @param ... Not in use.
#' @useDynLib rbvs mcplus_coef_gaussian_r mcplus_coef_binomial_r
#' @author Rafal Baranowski, Patrick Breheny
#' @export 
#' @references Zhang, Cun-Hui. "Nearly unbiased variable selection under minimax concave penalty." The Annals of Statistics (2010): 894-942.
#' @references Breheny, Patrick, and Jian Huang. "Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection." The Annals of Applied Statistics 5.1 (2011): 232.

mcplus.coef <- function(x,y,subsamples,
                              nonzero=NULL,
                              family=c("gaussian", "binomial"),
                              alpha=1.0,
                              gamma=3.0,
                              maxit=500,
                              tol=1e-2,
                              lambda.ratio=1e-6,
                              nlam=25,
                              ...){
  
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  y <- as.numeric(y)
  storage.mode(y) <- "double"
  
  subsamples <- as.matrix(subsamples)
  storage.mode(subsamples) <- "integer"
  
  
  #check if the sizes of x and y match
  
  n <- length(y)
  if(nrow(x) != n) stop("The number of rows in 'x' must be equal to the length of 'y'.")
  
  #check for NA's
  if(any(is.na(x))) stop("Matix 'x' cannot contain NA's.")
  if(any(is.na(y))) stop("Vector 'y' cannot contain NA's.")
  if(any(is.na(subsamples))) stop("Matrix 'subsamples' cannot contain NA's.")
  
  #check if subsamples are actually subsamples
  min.ind <- min(subsamples)
  max.ind <- max(subsamples)
  if(min.ind < 1 || max.ind > n) stop("Elements of 'subsamples' must be between 1 and length(y).")
  
  alpha <- as.double(alpha)[1]
  
  if(alpha < 0 || alpha > 1) stop("Parameter 'alpha' has to be between 0 and 1.")
  
  gamma <- as.double(gamma)[1]
  if(gamma < 1) stop("'gamma' must be greater than 1.")
  
  maxit <- as.integer(maxit)[1]
  
  if(maxit < 50){
    warning("Paramater 'maxit' must be greater than or equal to 50. Setting 'maxit' to 50.")
    maxit <- 50
  } 
  
  tol <- as.double(tol)[1]
  
  if(tol <= 0) stop("'tol' has to be positive.") 
  
  
  lambda.ratio <- as.double(lambda.ratio)[1]
  if(lambda.ratio < 0 || lambda.ratio > 1) stop("'lambda.ratio' has to be between 0 and 1.")
  
  nlam <- as.integer(nlam)[1]
  
  if(nlam < 2) {
    warning("Paramater 'nlam' must be greater than or equal to 50. Setting 'nlam' to 2.")
    maxit <- 2
  }
  
  
  #check if subsamples are actually subsamples
  min.ind <- min(subsamples)
  max.ind <- max(subsamples)
  
  
  
  m <- nrow(subsamples);
  B <- ncol(subsamples);
  p <- ncol(x)
  n <- nrow(x)
  
  if(is.null(nonzero)) nonzero <- m-1;
  
  family <- match.arg(family)
  
  if( family=="gaussian"){
    return(.Call("mcplus_coef_gaussian_r", 
                 as.matrix(subsamples),
                 as.matrix(x),
                 as.numeric(y),
                 as.integer(nonzero),
                 as.double(alpha),
                 as.double(gamma),
                 as.double(tol),
                 as.integer(maxit),
                 as.double(lambda.ratio),
                 as.integer(nlam)))   
  }else{
    return(.Call("mcplus_coef_binomial_r", 
                 as.matrix(subsamples),
                 as.matrix(x),
                 as.numeric(y),
                 as.integer(nonzero),
                 as.double(alpha),
                 as.double(gamma),
                 as.double(tol),
                 as.integer(maxit),
                 as.double(lambda.ratio),
                 as.integer(nlam)))
  }
  
}

#' @title Measure an impact of the covariates on the response using Pearson correlatio.
#' This function evaluates the Pearson correlation coefficient between the response \code{y} and each column in the design matrix \code{x} over subsamples in \code{subsamples}.
#' @param x Matrix with \code{n} observations of \code{p} covariates in each row.
#' @param y Response vector with \code{n} observations.
#' @param subsamples Matrix with \code{m} indices of \code{N} subsamples in each column. 
#' @param ... Not in use.
#' @return Numeric \code{p} by \code{N} matrix with Pearson correlations evaluated for each subsample.
#' @useDynLib rbvs pearson_cor_r
#' @export

pearson.cor <- function(x,y,subsamples,  ...){
  
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  y <- as.numeric(y)
  storage.mode(y) <- "double"
  
  subsamples <- as.matrix(subsamples)
  storage.mode(subsamples) <- "integer"

  
  #check if the sizes of x and y match
  
  n <- length(y)
  if(nrow(x) != n) stop("The number of rows in 'x' must be equal to the length of 'y'.")
  
  #check for NA's
  if(any(is.na(x))) stop("Matix 'x' cannot contain NA's.")
  if(any(is.na(y))) stop("Vector 'y' cannot contain NA's.")
  if(any(is.na(subsamples))) stop("Matrix 'subsamples' cannot contain NA's.")

  
  #check if subsamples are actually subsamples
  min.ind <- min(subsamples)
  max.ind <- max(subsamples)
  
  if(min.ind < 1 || max.ind > n) stop("Elements of 'subsamples' must be between 1 and length(y).")
  
  return(.Call("pearson_cor_r", 
               subsamples,
               x,
               y))  

} 

#' @title Measure an impact of the covariates on the response using the distance correlation 
#' This function evaluates the distance correlation between the response \code{y} and each column in the design matrix \code{x} over subsamples in \code{subsamples}.
#' @useDynLib rbvs distance_cor_r
#' @param x Matrix with \code{n} observations of \code{p} covariates in each row.
#' @param y Response vector with \code{n} observations.
#' @param subsamples Matrix with \code{m} indices of \code{N} subsamples in each column. 
#' @param index Positive scalar.
#' @param ... Not in use.
#' @return Numeric \code{p} by \code{N} matrix with distance correlations evaluated for each subsample.
#' @references   Maria L. Rizzo and Gabor J. Szekely (2014). energy: E-statistics
#' (energy statistics). R package version 1.6.1 (\url{http://CRAN.R-project.org/package=energy}).
#' @export

distance.cor <- function(x,y, subsamples, index=1, ...){
  
  
  #guarantee data types correct for C function
  
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  y <- as.numeric(y)
  storage.mode(y) <- "double"
  
  subsamples <- as.matrix(subsamples)
  storage.mode(subsamples) <- "integer"
  
  index <- as.numeric(index)
  
  #check if the sizes of x and y match
  
  n <- length(y)
  if(nrow(x) != n) stop("The number of rows in 'x' must be equal to the length of 'y'.")
  
  #check for NA's
  if(any(is.na(x))) stop("Matix 'x' cannot contain NA's.")
  if(any(is.na(y))) stop("Vector 'y' cannot contain NA's.")
  if(any(is.na(subsamples))) stop("Matrix 'subsamples' cannot contain NA's.")
  if(any(is.na(index))) stop("Scalar 'index' cannot contain NA's.")
  
  #check if subsamples are actually subsamples
  min.ind <- min(subsamples)
  max.ind <- max(subsamples)
  
  if(min.ind < 1 || max.ind > n) stop("Elements of 'subsamples' must be between 1 and length(y).")
  
  if(length(index) > 1){
    warning("Variable 'index' contains too many elements. Taking just the first one.")
    index <- index[1]
  }else if (length(index) <1) warning("Variable 'index' contains too few elements. Setting 'index' to 1.")
  
  if(index < 0) stop("Variable 'index' must be positive.")
  
  
  # Call C function calculating correlations

	return(.Call("distance_cor_r", 
	   subsamples,
	   x,
	   y,
     as.numeric(index)))
} 