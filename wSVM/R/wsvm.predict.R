wsvm.predict <- function(X, Y, new.X, new.Y, model, comp.error.rate = FALSE) {
# Description 
#     predict new.Y at new.X and compute error rate
# Usage
#     pred = wsvm.predict(X, Y, new.X, new.Y, model, comp.error.rate = F)
# Input
#     X = training X
#     Y = training Y
#     new.X = new X
#     new.Y = new Y
#     alpha = estimated coefficients (nrow(X) by 1)
#     bias = bias term
# Output
#     predicted.values = fitted value at new.X
#     g = sign(predicted.values)
#     error.rate = misclassification error rate

    Sign <- function(x){
       Sign <- sign(x)
       Sign[Sign == 0] <- 1
       return(Sign)
    }

   if(is.matrix(X)) X <- as.matrix(X)
   if(is.matrix(Y)) Y <- as.matrix(Y)
   if(is.matrix(new.X)) new.X <- as.matrix(new.X)
   
   K <- wsvm.kernel(X, new.X, model$kernel)
   predicted.values <- model$bias + t(Y * model$alpha) %*% K
   
   if (comp.error.rate == FALSE){
      predicted.Y = Sign(predicted.values)
      predict = list(predicted.values = predicted.values, predicted.Y = predicted.Y, error.rate = NULL)
   }
   else{
      predicted.Y = Sign(predicted.values)
      new.Y = as.vector(new.Y)
      error.rate = sum(predicted.Y != new.Y) / length(new.Y)
      predict = list(predicted.values = predicted.values, 
                     predicted.Y = predicted.Y, 
                     error.rate = error.rate)
   }        
   return(predict)   
}
