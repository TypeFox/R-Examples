#### Author: Adriano Zanin Zambom
#### contact: adriano.zambom@gmail.com
#### last modified: 2/Oct/2015
####
#### papers of reference: 
#### Jianqing Fan and Irene Gijbels (1995). Data-Driven Bandwidth Selection in Local Polynomial Fitting: Variable Bandwidth and Spatial Adaptation. Journal of the Royal Statistical Society. Series B (Methodological), Vol. 57, No. 2, pp. 371-394.
#### Masry (1996) MULTIVARIATE LOCAL POLYNOMIAL REGRESSION FOR TIME SERIES: UNIFORM STRONG CONSISTENCY AND RATES, J. Time Series Analysis, vol. 17 (November 1996), pp. 571-599.
#### Ruppert and Wand (1994) Multivariate locally weighted least squares regression, Annals of Statistics
 


##########################################################################################################################
# Function: localpoly.reg
#
##########################################################################################################################
localpoly.reg <- function(X, Y, points = NULL, bandwidth = "CV", gridsize = 30, degree.pol = 0, kernel.type = "epanech", deriv = 0) UseMethod("localpoly.reg")


print.localpoly.reg <- function(x,...)
{

   cat("Call:\n")
   print(x$call)    
   
    if (is.null(dim(x$x)))
    {
        plot(x$x,x$y, xlab="X",ylab="Y")
        lines(x$points[order(x$points)],x$predicted[order(x$points)])
    } 

    cat("\nbandwidth: ")
    cat(x$bandwidth)
    cat("\n\npredicted: ")
    cat(x$predicted)
    cat("\n")
}



##########################################################################################################################
# Function: localpoly.reg
# Parameters Received:
#
# X = Matrix with columns being the variables and rows the observations;
# Y = vector of observations;
# points = matrix with same number of columns of X, containing the (rows) points where the regression should be smoothed;
# bandwidth; degree.pol; kernel.type; deriv = which derivative, 0 is for the smoothed regression
#
# Values Returned:
# test = list
#   test$x = X = same as input
#   test$y = Y = same as input
#   test$points = points where the curve was estimated
#   test$bandwidth = used for the fit
#   test$predicted = predicted values
#
# Call: a function in C:local_poly_estimator
##########################################################################################################################
localpoly.reg.default <- function(X, Y, points = NULL, bandwidth = "CV", gridsize = 30, degree.pol = 0, kernel.type = "epanech", deriv = 0)    
{

    if ((!identical(kernel.type,"box")) && (!identical(kernel.type,"trun.normal")) && (!identical(kernel.type,"gaussian")) && (!identical(kernel.type,"epanech")) && (!identical(kernel.type,"biweight")) && (!identical(kernel.type,"triweight")) && (!identical(kernel.type,"triangular")))
       stop("\n\nInvalid kernel.type: ",kernel.type,"\n\n")
    if (!is.null(dim(Y)))
       stop("\n\n Y must be a univariate vector of observations \n\n")

    kernel.type = match(kernel.type,c("box","trun.normal","gaussian","epanech","biweight","triweight","triangular"))
    ## this will attribute a number to kernel.type corresponding to which string was inputed by the user: 
    ## "box" = 1, "trun.normal" = 2, etc

    if ((!is.numeric(gridsize)) || (gridsize < 2))
       stop("\n\nInvalid gridsize: gridsize must be an integer larger than 1.\n\n")
    
    if (!is.null(dim(X)) && (degree.pol > 2)) # X is a matrix
       stop("\n\npolynomial regression with degree > 2 available only for univariate X.\n\n") else
    if ((degree.pol < 0) || (!is.numeric(degree.pol)))
       stop("\n\ninvalid degree of polynomial.\n\n") else
    if ((deriv < 0) || (!is.numeric(deriv)))
       stop("\n\ninvalid derivative.\n\n")

   
    if (is.null(points)) 
       points = X else
    if (((is.null(dim(X))) && (!is.null(dim(points)))))
       stop("\n\nX has dimension 1, points cannot be of dimension ",dim(points),", procedure stopped.\n\n") else
    if ((is.null(dim(points))) && ((!is.null(dim(X))))) # X is a matrix and points is a vector
    {
       if (dim(X)[2] != length(points))
          stop("\n\n 'points' is required to be a matrix or vector with same dimension of X, with ",dim(X)[2]," covariates, procedure stopped.\n\n")
    } else
    if ((!is.null(dim(points))) && ((!is.null(dim(X))))) # X and points are matrices
    {
       if (dim(X)[2] != dim(points)[2])
          stop("\n\n 'points' is required to be a matrix or vector with same dimension of X, with ",dim(X)[2]," covariates, procedure stopped.\n\n") 
    }

   
    if (is.numeric(bandwidth)) ## bandwidth is a number of a vector of numbers
    {
       if (is.null(dim(X)))  ## X is a vector -> bandwidth must be vector or single number
       {
          if (length(bandwidth) == 1) ## single number
          {
             if (bandwidth <= 0)
                stop("\n\n Argument bandwidth must be positive. \n\n")
          } else ## vector
          if ((length(bandwidth) != length(points)) || (is.matrix(bandwidth)) || (min(bandwidth) <= 0))
             stop("\n\n An input bandwidth must be a number or positive vector with the same length as 'points'. \n\n")
       } else ## X is a matrix
       if ((length(bandwidth) != length(points)) && (dim(X)[2] != length(bandwidth)))
          stop("\n\n For multivariate X, bandwidth must be a vector of the same dimension of the columns of X, columns of 'points' or rows of 'points'. \n\n")
    } else # bandwidth is a string
    if ((!identical(bandwidth,"CV")) && (!identical(bandwidth,"GCV")) && (!identical(bandwidth,"CV2")) && (!identical(bandwidth,"GCV2")) && (!identical(bandwidth,"Adp")))
       stop("\n\n Invalid bandwidth: ",bandwidth,".\n\n") else
    {
        if ((identical(bandwidth,"Adp")) && (!is.null(dim(X))))
           stop("\n\n Adaptive bandwidth choice is only valid for univariate X.\n\n")
    
        if (is.null(dim(X))) # X is a vector
        {
           if (identical(bandwidth,"CV") || identical(bandwidth,"CV2"))
              bandwidth = 0 else
           if (identical(bandwidth,"GCV") || identical(bandwidth,"GCV2"))
              bandwidth = -1 else
              bandwidth = -4
        } else # X is a matrix
        {
            if (identical(bandwidth,"CV"))
               bandwidth = 0 else
            if (identical(bandwidth,"GCV"))
               bandwidth = -1 else
            if (identical(bandwidth,"CV2"))
               bandwidth = -2 else
            if (identical(bandwidth,"GCV2"))
               bandwidth = -3
        }
        
    }

    
    ## create matrix of bandwidths:
    ## repeated for each data point of 'points' if matrix 
    ## of bandwidths has not been not specified by user
    banda = bandwidth
    if (bandwidth[1] > 0) # bandwidth is a number of a vector of values
    {
       if (!is.null(dim(X))) # X is a matrix
       {
          if (length(points) != dim(X)[2]) ## points is a matrix
          if (length(bandwidth) != length(points))
          {
             banda = matrix(0,dim(points)[1],dim(points)[2])
             for (i in 1:dim(points)[2])
                banda[,i] = rep(bandwidth[i],dim(points)[1])
          }
       } else # X is a vector
         if (length(bandwidth) == 1)
            banda = rep(bandwidth,length(points))
    } else # bandwidth is 0, -1, -2, -3, -4
    {
        if (is.null(dim(points)))
            banda = rep(bandwidth[1],length(points)) else
            banda = matrix(bandwidth[1],dim(points)[1],dim(points)[2])
    }
    


    
    
    ## construct grid for cross-validation
    grid = 1  # wont be used
    if ((bandwidth[1] == 0) || (bandwidth[1] == -1) || (bandwidth[1] == -2) || (bandwidth[1] == -3))
    {    
       if (is.null(dim(X)))
          grid = seq(quantile(dist(X),.01),max(dist(X))/2,length=gridsize) else
       {
           grid = matrix(0,gridsize,dim(X)[2])
           for (i in 1:dim(X)[2])
              grid[,i] = seq(quantile(dist(X[,i]),.01),max(dist(X[,i]))/2,length=gridsize)
       }
    }
      
      
    if (bandwidth[1] == -4) # adaptive bandwidth, binned
    {
           ordem = order(X) # need to order X to find correct bins
           X = X[ordem]
           Y = Y[ordem]
           points = points[ordem]
           bins = floor(1.5*length(X)/(10*log(length(X))))
           n_bins = floor(length(X)/bins)
           
           band = 0
           for (i in 1:bins)  ## find band with cross-validation in each bin
           {
               grid = seq(quantile(dist(X[(n_bins*(i-1)+1):(n_bins*i)]),.01),max(dist(X[(n_bins*(i-1)+1):(n_bins*i)]))/2,length=gridsize)
               m1 = .Call("local_poly_estimator",t(X[(n_bins*(i-1)+1):(n_bins*i)]),Y[(n_bins*(i-1)+1):(n_bins*i)], t(X[(n_bins*(i-1)+1):(n_bins*i)]), t(rep(0,length=length(X))), t(grid), degree.pol, kernel.type, deriv) 

               band[i] = m1$bandwidth[1]
           }
           banda = rep(band[1],n_bins)
           for (i in 2:bins)
              banda = c(banda, rep(band[i],n_bins))
           if ((bins*n_bins) != length(X))
              banda = c(banda, rep(band[bins],(length(X) - (n_bins*bins))))
           
           grid = seq(1,length(X)/2,length=gridsize)
           m1 = .Call("local_poly_estimator",t(seq(1,length(X))), banda, t(seq(1,length(X))), t(rep(0,length=length(X))), t(grid), degree.pol, kernel.type, deriv) 
           banda = m1$predicted
    }

    # in C, a n by d matrix is a vector of size (nd)
    # I send the transpose so that in C, the first d elements are the 1rst row, and so on
    m1 = .Call("local_poly_estimator",t(X),Y, t(points), t(banda), t(grid), degree.pol, kernel.type, deriv) 
                
    
   test = list()  ## test is the output list of values of the function
   test$x = X
   test$y = Y
   test$points = points
    
    
   if (is.null(dim(X)))
   {
      if ((length(bandwidth) == length(points)) || (bandwidth[1] == -4))
         test$bandwidth = m1$bandwidth
      else
         test$bandwidth = m1$bandwidth[1]
   } else
   if ((bandwidth[1] == 0) || (bandwidth[1] == -1) || (bandwidth[1] == -2) || (bandwidth[1] == -3) || (length(bandwidth) == dim(X)[2]))
        test$bandwidth = m1$bandwidth[1:dim(X)[2]]
     else
        test$bandwidth = m1$bandwidth
    
    
    
    
   test$predicted = m1$predicted

   result = test
   result$call <- match.call()

   class(result) <- "localpoly.reg"
   result

}





