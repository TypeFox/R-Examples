`mu.varknown` <-
function(len,lambda,n0,level=0.95)
{
  # Returns the minimal sample size to ensure an average prob. coverage of
  # 'level' with a symmetric region about the posterior mean
  # when sampling from a normal distribution with known variance (Adcock, 1988)

  # Using the same notation as Adcock
  var <- 1/lambda
  d <- len/2

  m <- n0
  m <- ceiling(((qnorm((level+1)/2)/d)^2)*var-m)
  
  max(0, m)
}

