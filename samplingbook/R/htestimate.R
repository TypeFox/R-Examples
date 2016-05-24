htestimate <-
function(y, N, PI, pk, pik, method="yg"){

# input treatment
    if(missing(y)) stop("Wrong input: ", sQuote("y")," has to be given as a vector of observations")
    else n <- length(y)
    if(!(method == 'yg' || method == 'ht' || method == 'hh' || method == 'ha')){
      stop("Wrong input: Only types ", sQuote("yg")," (Yates and Grundy), ", sQuote("ht")," (Horvitz-Thompson), ", sQuote("ha")," (Hajek) or ", sQuote("hh")," (Hansen-Hurwitz).")
    }
    if(missing(N) || !is.numeric(N)) stop("Wrong input: Population size ", sQuote("N")," has to be a positive integer smaller than ", sQuote("Inf"),".") 
    if(method == 'yg' || method == 'ht'){
      if(missing(PI) || !is.matrix(PI)) stop("Wrong input: For using methods ", sQuote("yg")," or ", sQuote("ht")," matrix ", sQuote("PI")," has to be given.")
      if(dim(PI)[1]!=dim(PI)[2]) stop("Wrong input: dimension error, ", sQuote("PI")," has to be a square matrix.") 
      if(dim(PI)[1] != length(y)) stop("Wrong input: number of rows or columns of ", sQuote("PI")," has to be as sampling size ", sQuote("n"),".")
      pk <-  diag(PI)
    }  
    if(method == 'hh' || method == 'ha'){
      if(missing(pk) || !is.vector(pk)) stop("Wrong input: For using methods ", sQuote("hh")," or ", sQuote("ha")," vector ", sQuote("pk")," of inclusion probabilities has to be given.")   
      if(length(pk) != n) stop("Wrong input: vector ", sQuote("pk")," of inclusion probabilities has to have length of sample size ", sQuote("n"),".")
      if(method == 'ha'){
        if(missing(pik)) warning("Without input of ", sQuote("pik")," just approximative calculation of Hajek method is possible.")
        else{
          if(length(pik) != N ) stop("Wrong input: vector ", sQuote("pik")," of inclusion probabilities has to have length of population size ", sQuote("N"),".")
          if(!all(pk %in% pik) ) stop("Wrong input: vector ", sQuote("pk")," of inclusion probabilities has to be subset of ", sQuote("pik"),".")
        }
      }                 
    } 
    if(missing(pik)) pik <- NULL
# mean calculation
  mean.estimate <-  sum(y/pk)/N
# variance calculations
  if(method=="yg" || method=="ha"){
    if (method == "ha"){                                 
      if(is.null(pik)) dpik <- sum((1-pk)) # approximative: length(pk)=n
      else dpik <-  sum(pik*(1-pik)) # exact: length(pik)=N
      PI <- pk %*%  t(pk)*( 1- ( (1-pk) %*% t(1 - pk)/ dpik))
      diag(PI) <- pk
    }                                          
    A <- ((matrix(pk) %*% t(matrix(pk)) - PI)/PI)*(kronecker(y/pk,matrix(1,1,n)) - kronecker(matrix(1,n,1),t(y/pk)))^2
    var.estimate <- 1/(2*N^2)*sum(A)
  }
  if (method=="ht"){
    A <- ( PI -  ( pk %*% t(pk)))/( PI * ( pk %*% t(pk))) * ( y %*% t(y))
    diag(A) <-  (1 - pk)/pk^2 * y^2
    var.estimate <- 1/(N^2) * sum(A)
  }
  if (method == "hh"){                                                           
    pkn <-  pk/n
    A <- (y/(N*pkn) - mean.estimate)^2
    var.estimate <- sum(A)/(n*(n-1))
    PI <- NULL
  }
  
# se calculation
  if(var.estimate>0) se <- sqrt(var.estimate)
  else{
    se <- NA
    warning("Standard error is ", sQuote("NA"),", because calculations for variance of mean has been not positive.")
  }  
    
# return argument  
  ret <- list()
  ret$call <- list(y=y,N=N,PI=PI,pk=pk,pik=pik,method=method)
  ret$mean <- mean.estimate
  ret$se <- se
  structure(ret, class="htestimate")
}
