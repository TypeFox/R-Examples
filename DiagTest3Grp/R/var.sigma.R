var.sigma <-
function(s,n)
  {
   ###This function calculates the variance of the MLE estimator on a normal SD , i.e., (sample SD)
   ### var(s)=var(sqrt(sum(x_i-x_bar)^2)/(n-1))
   ###s: population SD
   ##n: sample size
    0.5*s^2/(n-1)
  }

