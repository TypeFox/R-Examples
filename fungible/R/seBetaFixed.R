#------------------------------------------------------#
# This function computes the Cov of std Reg coeffs for #
# fixed X  using a method described by Yuan and Chan   #
# (2011)                                               #
#                                                      #
#                                                      #
# Arguments:                                           # 
# X            - matrix of predictor scores            #  
# y            - vector of criterion scores            #   
# cov.x        - covariance matrix for predictors      #    
# cov.xy       - vector of covariances between         #  
#                predictors and criterion              #  
# var.y        - criterion variance                    #  
# var.error    - var.error                             #
# Nobs         - number of observations                # 
#                                                      #
# This function accepts either (1) raw data, or (2)    #
# second-order moments (covariances or correlations)   #
# and sample size.                                     #
#                                                      #
# Output                                               #
#  seBeta      - normal theory standard errors for     #
#                standarized regression coefficients   #
#                with Fixed predictors.                #
#  covBeta     - normal theory covariance matrix of    #
#                standarized regsion coefficients      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



seBetaFixed <- function(X=NULL,y=NULL,cov.x=NULL,cov.xy=NULL,
	                    var.y=NULL,var.error = NULL,Nobs=NULL) {

########################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~ Error Checking ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  if(is.null(X) & !is.null(y)) 
    stop("\n y is not defined\n Need to specify both X and y\n")
  
  if(!is.null(X) & is.null(y)) 
    stop("\n X is not defined\n Need to specify both X and y\n")
			
  if(is.null(X) & 	is.null(y)) {
	
    if(is.null(cov.x) | is.null(cov.xy) | is.null(var.y) | is.null(Nobs))
      stop("\nYou need to specify all covariances and the number of observations\n")
      
    N <- Nobs	
    p <- nrow(cov.x)
		
  } else {
    
# if function is given data, compute all the covariances and variances    
	
    cov.x <- cov(X)
    cov.xy <- cov(X,y)
    var.y <- var(y)
    N <- length(y)
    p <- ncol(X)

  }

# create diagonal matrix of predictor standard deviations (Dx),
# unstandardized regression coefficients (b), and
# error variance from regression model (var.error)

  var.y <- as.numeric(var.y)
  Dx <- diag(sqrt(diag(cov.x)))	

  b <- solve(cov.x)%*%cov.xy
  
# if error variance is not supplied, calculate it
  
  if(is.null(var.error)) var.error <- as.numeric(var.y - t(b)%*%cov.x%*%b)*(N-1)/(N-p-1)	

# create Gamma matrix which is the 
# covariance matrix of b (unstandardized) and sigma^2 y (variance of y)

  var.b <- solve(cov.x)*var.error
  cov.b.var.y <- 2*b*var.error  
  var.var.y <- 4*(t(b)%*%cov.x%*%b)*var.error + 2*var.error^2
  Gamma <- rbind(cbind(var.b,cov.b.var.y),c(cov.b.var.y,var.var.y))

# Create jacobian for the function g(b,var.y) = Dx%*%b/sd.y

  deriv_b <- Dx*(1/sqrt(var.y))
  deriv_var.y <- -Dx%*%b/(2*var.y^(3/2))	
	
  jacobian <- cbind(deriv_b, deriv_var.y)	
	
# create covariance matrix of standardized regression coefficients
  
  cov.mat <- jacobian%*%Gamma%*%t(jacobian)/N

# name the covariance matrix  
  
  b.nms <- NULL
  for(i in 1:p) b.nms[i] <- paste("b",i,sep='')

  rownames(cov.mat) <- colnames(cov.mat) <- b.nms
  
  list(covBeta=cov.mat, seBeta=sqrt(diag(cov.mat)))	
	
}	

