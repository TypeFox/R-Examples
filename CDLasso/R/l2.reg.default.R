l2.reg.default <-
function (X, Y, lambda=1)
{
  
 # if ((!is.numeric(lambda))|lambda<=0)
 #   stop ("lambda should be a positive number")
  
  if (length(Y)!=ncol(X))
    stop("Dimention doesn't match! 
         Columns of feature matrix X must be the number of cases")
  
  # initialize input
  cases <- length(Y)
  predictors <- nrow(X)
  # add intercept row
  X <-rbind(rep(1,ncol(X)),X)
  
  return_data <- .Fortran ("LASSO_PENALIZED_L2_REGRESSION", 
                           as.double(X),
                           Y=as.double(as.vector(Y)),
                           as.double (lambda),
                           as.integer (cases),
                           as.integer (predictors+1),
                           L2 = as.double(0),
                           r = as.double(as.vector(rep(0,cases))),
                           objective = as.double(0),                        
                           penalty = as.double(0),
                           estimate =  as.double(as.vector(rep(0,predictors+1))),
                           PACKAGE = "CDLasso")
  
  selected<-c()
  nonzeros<-0
  for (j in 2:(predictors+1))
  {
    if (return_data$estimate[j]!=0) 
    {
      nonzeros<-nonzeros+1
      selected<-c(selected,j-1)
    }
  }
  
  out <- list (X = X,
  			  Y = Y,
  			  intercept = return_data$estimate[1],
              estimate = return_data$estimate[2:(predictors+1)],	               
              cases = cases,
              predictors = predictors,
              lambda = lambda,
              L2 = return_data$L2,
              residual = return_data$r,
              nonzeros = nonzeros,
              selected = selected,
              call=sys.call ())
  class (out) <- "l2.reg"
  
  return(out)
}
