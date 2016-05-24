lqa.control <-
function (x = NULL, var.eps = .Machine$double.eps, max.steps = 5000, conv.eps = 1e-03, conv.stop = TRUE, c1 = 1e-08, digits = 5, ...)
{
###  x = NULL  ... object of class 'lqa'. This optional argument is just included to be in line with the S3 class concept
###  var.eps   ... tolerance in checking for zero variance of some regressors
###  max.steps ... maximum number of steps in the P-IRLS respective Boosting algorithms
###  conv.eps  ... tolerance for convergence break in parameter updating
###  conv.stop ... whether or not to stop the iterations when estimated coefficients are converged
###  digits    ... number of digits to round the tuning parameter candidates as names in the 'loss.array' (in function 'cv.lqa')
###  ...       .... further arguments 

   if (!is.numeric (var.eps) || var.eps < 0)
     stop ("value of var.eps must be >= 0")

   max.steps <- as.integer (max.steps)
   digits <- as.integer (digits)

   if (max.steps < 0)
     stop ("max.steps must be positive integer")

   if (!is.numeric (conv.eps) || conv.eps <= 0)
     stop ("value of conv.eps must be > 0")
  
   if (!is.logical (conv.stop))
     stop ("conv.stop must be 'TRUE' or 'FALSE'")

   if (!is.numeric (c1) || c1 < 0)
     stop ("value of 'c1' must be >= 0")

   if (!is.numeric (digits) || digits < 0)
     stop ("value of 'digits' must be >= 0")

   list (var.eps = var.eps, max.steps = max.steps, conv.eps = conv.eps, conv.stop = conv.stop, digits = digits, c1 = c1)
}

