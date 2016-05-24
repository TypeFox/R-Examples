#### Author: Adriano Zanin Zambom
#### contact: adriano.zambom@gmail.com
#### last modified: 2/Oct/2015
####
#### papers of reference: 
#### Zambom, A. Z. and Akritas, M. G. (2012). a) Nonparametric Model Checking and Variable Selection. Statistica Sinica, v. 24, pp. 1837.
#### Zambom, A. Z. and Akritas, M. G. (2012). b) Signicance Testing and Group Variable Selection. Journal of Multivariate Analysis, v. 133, pp. 51.




#--------------------------------------------------------------------------------
# this function is used when we dont need to estimate the residuals from a fit,
# but instead, just use the ANOVA test for X univariate
#--------------------------------------------------------------------------------
anova_test_univariate <- function(X, Y, p)   # include the edges
{
   n = length(X)
   xi = Y
   xi = xi[order(X)]

   V = matrix(0, n-p+1, p)
   for(ind in ((p+1)/2):(n-(p+1)/2+1))                # I am excluding the edges
      V[ind-(p+1)/2+1,] = xi[(ind-(p+1)/2+1):(ind+(p+1)/2-1)]

   nv = dim(V)[1]
   t_mean = mean(V)

   MST = sum((rowMeans(V) - t_mean)^2)*p/(nv-1)
   MSE = sum((V - rowMeans(V))^2)/(nv*p - nv)

   tau = sum((xi[2:(n-2)] - xi[1:(n-3)])^2*(xi[4:n] - xi[3:(n-1)])^2)/(4*(n-3))
   if (tau == 0) ### when the estimate interpolates, tau = 0, and an error would occur, so we set it > 0
      tau = 0.0000000000001

   p_value = (1-pnorm((sqrt(nv)*(MST-MSE))/(sqrt(2*p*(2*p-1)*tau/(3*(p-1))))))

   return(p_value)
}


##########################################################################################################################
# Function: npmodelcheck
#
##########################################################################################################################

npmodelcheck <- function(X, Y, ind_test, p = 7, degree.pol = 0, kernel.type = "epanech", bandwidth = "CV", gridsize = 30, dim.red = c(1,10)) UseMethod("npmodelcheck")


print.npmodelcheck <- function(x,...)
{
    cat("Call:\n")
    print(x$call)
    cat("\nP-value of the test: ")
    options(digits=22)
    cat(x$p_value, "\n")
}

summary.npmodelcheck <- function(object, ...)
{

    res <- list(call=object$call,
                b = object$bandwidth, p = object$p_value, pred = object$predicted)
    class(res) <- "summary.npmodelcheck"
    res 
}


print.summary.npmodelcheck <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    options(digits=22)
    cat("\nbandwidth: ")
    cat(x$b)
    cat("\n\nP-value of the test: ")
    cat(x$p)
    cat("\n\nPredicted: ")
    cat(x$pred)
    cat("\n")
}




##########################################################################################################################
# Function: npmodelcheck.default
# Parameters Received:
#
# X = Matrix with columns being the covariates and rows the observations
# Y = vector of observations
# ind_test = index of column of X to be tested 
# p = number of covariate values included in the window W_i, has to be an odd number!
# degree.pol = the degree of the polynomial to fit the nonparametric regression
# kernel.type = type of kernel (box, truncated normal, etc) for estimation of m_1
# bandwidth = "CV", "GCV", "CV2", "GCV2", "Adp", number, vector or matrix
# dim.red: vector where 1st entry indicates 1 for SIR and 2 for SPC. Second entry indicates number of slices (if SIR) or number of components (if SPC)
#
# Values Returned:
# test = list containing the bandwidth used, the predicted Y values and the p-value of the test 
##########################################################################################################################
npmodelcheck.default <- function(X, Y, ind_test, p = 7, degree.pol = 0, kernel.type = "epanech", bandwidth = "CV", gridsize = 30, dim.red = c(1,10))
{
   
   if (!is.null(dim(Y))) stop(paste("Y must be a vector"))
   if ((p%%2 == 0) || (p < 2) || (p >= length(Y)) || (!is.numeric(p))) stop(paste("p must be an odd number greater or equal to 3 and smaller than the number of observations."))
   if (is.null(dim(X)))
   { 
      if (ind_test != 1) stop(paste("enter a valid covariate to be tested")) 
      if (length(X) != length(Y)) stop(paste("length of covariate vector and response must be the same")) 
      test = list(bandwidth = NULL, predicted = NULL, p_value = anova_test_univariate(X, Y, p))
   }  else   
   {
   if (dim(X)[1] < dim(X)[2]) stop(paste("Number of variables larger than observations")) else
   if (length(X[,1]) != length(Y)) stop(paste("length of covariate vector and response must be the same")) else
   if (length(ind_test) == 1)    
      if (ind_test <= 0) stop(paste("enter a valid covariate number to be tested")) else  
      if (ind_test > dim(X)[2]) stop(paste("enter a valid covariate number to be tested")) 
   if (sum(as.numeric(ind_test < (dim(X)[2])+1)) < length(ind_test)) stop(paste("enter valid covariate indices to be tested"))  else
   if (sum(as.numeric(ind_test > 0)) < length(ind_test)) stop(paste("enter valid covariate indices to be tested"))  else
   if (length(ind_test) > (dim(X)[2])) stop(paste("Index of covariates to be tested must be less than the total number of covariates")) else   
   {
      n = dim(X)[1]
                                                         ## X_I is the matrix that will be passed to estimate m_1(X_I)
      X_I = X[,c(-ind_test,-which(diag(var(X)) == 0))]   ## it eliminates covariates that are constants (var = 0)

      if (dim.red[1] == 0)   ### No Dimension Reduction
      {
         X_I = X[,-ind_test]
      } else
      if (!is.null(dim(X_I)))
      if (dim(X_I)[2] > 1)
      if (dim.red[1] == 1)   ### SIR
      { 
         if (length(dim.red) != 2) stop(paste("incorrect input for dim.red! must be a vector with 2 entries")) else  
         if (dim.red[2] < 2) stop(paste("Number of slices must be greater than 1")) else  
         if (dim.red[2] > (n-1)) stop(paste("Number of slices must be less than the number of observations")) else  
         {
            slices = dim.red[2]
            s0 <- dr(Y ~ X_I, nslices = slices, method = "sir")

            chi_test = 0
            for (K in 0:min(slices-2,length(s0[10]$evalues)))
                chi_test[K+1] = 1-pchisq(n*(length(s0[10]$evalues)-K)*mean((sort(s0[10]$evalues))[1:(length(s0[10]$evalues)-K)]),(length(s0[10]$evalues)-K)*(slices-K-1))
            number_factors = length(which(chi_test < 0.05))
  
            if (number_factors > 1) 
               X_I =  t(t(s0$evectors[,1:number_factors])%*%t(X_I)) else
               X_I =  t(t(s0$evectors[,1:number_factors])%*%t(X_I))[,1] 
         }
      } else
      if (dim.red[1] == 2)   ### SPC
      {
         if (length(dim.red) != 2) stop(paste("incorrect input for dim.red! must be a vector with 2 entries")) else  
         if (dim.red[2] > dim(X_I)[2]) stop(paste("Number of Principal Components must be less than the dimension of X[-ind_test]")) else
         if (dim.red[2] <= 0) stop(paste("Number of Principal Components must be greater than 0")) else
         {
            X_I_aux = X_I
            s0 = 0
            for (prin in 1:(dim(X_I_aux)[2]))
               s0[prin] = anova_test_univariate(X_I_aux[,prin], Y, p)

            quais = which(s0 < 0.3)       #### theta = 0.3 threshold parameter for SPC
            if (length(quais) == 0) quais = 1    ### if no p-value is < theta, get the smallest
            quais = order(s0)[1:(length(quais))]   
            Princ_Comp = princomp(X_I_aux[,quais])
            if (length(Princ_Comp$loadings[,1:min(dim.red[2],length(quais))]) == 1)
               X_I =  as.numeric(X_I_aux[,quais]*(Princ_Comp$loadings[,1:min(dim.red[2],length(quais))])) else
               if (min(dim.red[2],length(quais)) == 1)
                  X_I =  as.numeric(X_I_aux[,quais]%*%(Princ_Comp$loadings[,1:min(dim.red[2],length(quais))])) else
                  X_I =  X_I_aux[,quais]%*%(Princ_Comp$loadings[,1:min(dim.red[2],length(quais))]) 
         }
          
      } else stop(paste("first entry of dim.red must be 1 = SIR or 2 = SPC")) 


      if (length(X_I)/n == 0)  ## this is when X[-ind_test] = NULL, that is, test all covariates of X
      {     
         m1 = list()
         m1$bandwidth = 0
         m1$predicted = mean(Y) 
         xi = Y - m1$predicted
      } else
      {
         m1 = localpoly.reg(X_I, Y, X_I, degree.pol = degree.pol, kernel.type = kernel.type, bandwidth = bandwidth, gridsize = gridsize)   
         xi = Y - m1$predicted
      }




      if (length(ind_test) == 1)      ### Z univariate
         Z = X[,ind_test] else
      {                               ### Z multivariate, use 1sp SPC
         s0 = 0
         for (prin in 1:length(ind_test))
            s0[prin] = anova_test_univariate(xi, X[,ind_test[prin]], p)
         quais = which(s0 < 0.3)       #### theta = 0.3 threshold parameter for SPC
         if (length(quais) == 0) quais = 1    ### if no p-value is < theta, get the smallest
         quais = order(s0)[1:(length(quais))]   
         Princ_Comp = princomp(X[,ind_test[quais]])
         if (length(quais) == 1)
            Z = as.numeric(X[,ind_test[quais]]%*%t(Princ_Comp$loadings[,1])) else
            Z = as.numeric(X[,ind_test[quais]]%*%(Princ_Comp$loadings[,1])) 
      }





      xi = xi[order(Z)]   ### organize windows by the covariate being tested
      V = matrix(0, n-p+1, p)
      for(ind in ((p+1)/2):(n-(p+1)/2+1))                # Construct V, excluding the edges
         V[ind-(p+1)/2+1,] = xi[(ind-(p+1)/2+1):(ind+(p+1)/2-1)]

      nv = dim(V)[1]
      t_mean = mean(V)

      MST = sum((rowMeans(V) - t_mean)^2)*p/(nv-1)
      MSE = sum((V - rowMeans(V))^2)/(nv*p - nv)

      tau = sum((xi[2:(n-2)] - xi[1:(n-3)])^2*(xi[4:n] - xi[3:(n-1)])^2)/(4*(n-3))
      if (tau == 0)
         tau = 0.0000000000001

      test = list()
      test$bandwidth = m1$bandwidth
      test$predicted = m1$predicted
      test$p_value = (1-pnorm((sqrt(nv)*(MST-MSE))/(sqrt(2*p*(2*p-1)*tau/(3*(p-1))))))
   }
   }


   result = test
   result$call <- match.call()

   class(result) <- "npmodelcheck"
   result
}


