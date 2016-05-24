logit.reg.default <-
function (X, Y, lambda=1)
{  
  
    
  if (is.vector(X)){X = t(as.matrix(X))}
  
  if (length(Y)!=ncol(X))
   	stop("Dimension doesn't match! 
         Columns of feature matrix X must be the number of cases")


         
  people <- length(Y)
  parameters <- nrow(X) + 1
  # add intercept row
  X <-rbind(1,X)
  X1<-t(X)
  
  return_data <- .Fortran ("LASSO_PENALIZED_ESTIMATION", 
                           as.double(as.matrix(X1,people,parameters)),
                           Y=as.double(as.vector(Y)),
                           as.double (lambda),
                           as.integer (people),
                           as.integer (parameters),
                           objective = as.double (0),
                           loglikelihood = as.double(0),
                           r = as.double(as.vector(rep(0,people))),                       
                           estimate =  as.double(as.vector(rep(0,parameters))),
                           PACKAGE = "CDLasso")
  
  selected<-c()
  nonzeros<-0
  for (j in 2:(parameters))
  {
    if (return_data$estimate[j]!=0)  
    {
      nonzeros<-nonzeros+1
      selected<-c(selected,j-1)
    }
  }
  
  out <- list (X = X,
  			  Y = Y,
  			  cases = people,
              predictors = parameters-1,
              lambda = lambda,
              probs = return_data$probs,
              objective = return_data$objective,
              loglikelihood = return_data$loglikelihood,
              residual = return_data$r,
              intercept = return_data$estimate[1],
              estimate = return_data$estimate[2:length(return_data$estimate)],
              nonzeros = nonzeros,
              selected = selected,
              call=sys.call ())
  class (out) <- "logit.reg"
  
  return(out)
}
