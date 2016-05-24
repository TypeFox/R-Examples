var.mu <-
function(s,n)
  {###The function calculates the variance of MLE estimator on sample mean mu_hat, i.e., (sample mean)
   ## E(mu_hat)=E(x_bar)=s^2/n
   ###s: population SD or its estimator, the samples drawn from a normal distribution 
   ##n: sample size
    s^2/n
  }

