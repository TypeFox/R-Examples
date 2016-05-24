lambda.check <-
function (lambda)
{
  if (is.null (lambda))  ## check for existence of attribute 'lambda'
    stop ("Tuning parameter lambda missing \n")
 
  if (any (lambda < 0))   ## check for non-negativity of lambda
    stop ("lambda must be non-negative \n")
}

