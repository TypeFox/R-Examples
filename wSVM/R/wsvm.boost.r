wsvm.boost <- function(X, Y, new.X, new.Y, c.n, B = 50, kernel.type = list(type = "rbf", par= 0.5), C = 4, eps = 1e-10, plotting = FALSE){

# Description 
#     Boost wsvm algorithm
# Usage
#     boo = wsvm.boost(X, Y, new.X, new.Y, c.n, B = 50, kernel.type = list(type = "rbf", par= 0.5), C = 4, eps = 1e-10, plotting = FALSE)
# Input
#     X = input variable matrix
#     Y = output variable vector which will be declared as a matrix in SVM
#     new.X = new X
#     new.Y = new Y
#     c.n = parameter of weight term
#     B = number of iteration
#     kernel.type = a type of kenrel
#     C = regulazation parameter
#     eps = epsilon
#     plotting = logical value for plotting
# Output
#     error.rate = error rates
#     predicted.model = predicted model

   if(!is.matrix(X)) X <- as.matrix(X)
   if(!is.matrix(Y)) Y <- as.matrix(Y)
   if(!is.matrix(new.X)) new.X <- as.matrix(new.X)
   if(!is.matrix(new.Y)) new.Y <- as.matrix(new.Y)

   n.data <- nrow(X)
   
   misc.rate <- rep(0, B)
   alpha <- rep(0, B)
   B.error.cumul <- rep(0, B)
   final.test <- rep(0, length(new.Y))
   
   for (i in 1 : B){ 
      W.fit <- wsvm(X, Y, c.n, kernel = kernel.type, C=C, eps = eps)  
      fit.S <- wsvm.predict(X, Y, X, Y, W.fit, comp.error.rate = TRUE) 
      fit.S.Y <- t(fit.S$predicted.Y)

      misc.rate[i] <- sum(c.n * (Y != fit.S.Y)) / sum(c.n)
      if (misc.rate[i] <= eps) break
      alpha[i] <- log((1 - misc.rate[i]) / misc.rate[i])
      
      ######## the higher misclassfication rate is , the lower log-value is 
      c.n <- c.n * exp(alpha[i] * (Y != fit.S.Y))
      
      test.W <- wsvm.predict(X, Y, new.X, new.Y, W.fit, comp.error.rate = TRUE)
      test.W.Y <- test.W$predicted.values

      final.test <- final.test + alpha[i] * test.W.Y

      if (misc.rate[i] >= 0.45)  c.n <- rep(1 / n.data , n.data)      

      last.test <- sign(final.test)
      last.test <- t(last.test)

      B.error <- Error.rate(last.test, new.Y)
      B.error.cumul[i] <- B.error

   }
   
   res <- list(error.rate = B.error.cumul, predicted.model = test.W)
   if(plotting) plot(B.error.cumul, ylim = c(0.2,0.35), type='l')
   return(res)
}
