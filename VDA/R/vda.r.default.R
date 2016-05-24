vda.r.default <-
function (x, y, lambda=1/length(y))
{
  
  if ((!is.numeric(lambda))|lambda<=0)
    stop ("lambda should be a positive number")
  
  if (length(y)!=nrow(x))
    stop("Dimention doesn't match! 
         Rows of feature matrix X must be the number of cases")
  
  # initialize input
  cases <- length(y)
  classes <- length(unique(y))
  features <- ncol(x)
  # add intercept col
  feature_i <- as.matrix(cbind(rep(1,nrow(x)),x))
#  colnames(feature_i)[1] <- "intercept"

  return_data <- .Fortran ("VDA", 
                           stand.feature = as.double (feature_i),
                           as.integer (as.vector (y)),
                           as.integer (cases),
                           as.integer (classes),
                           as.integer (features),
                           as.double (lambda),
                           predicted = as.integer (rep (0,cases)),
                           coefficient = as.double (matrix (0,classes-1,features+1)),
                           training_error_rate = as.double (0),
                           PACKAGE = "VDA")
  
  out <- list (feature = feature_i,
              stand.feature = matrix (return_data$stand.feature,cases,features+1),
              class = y,
              cases = cases,
              classes = classes,
              features = features,
              lambda = lambda,
              predicted = return_data$predicted,
              coefficient = matrix (return_data$coefficient,classes-1,features+1),
              training_error_rate = return_data$training_error_rate,
              call=sys.call ())
  class (out) <- "vda.r"
  
  return(out)
}
