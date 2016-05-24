`mu.acc` <-
function(len, alpha, beta, n0, level=0.95)
{
  # First make sure that prior knowledge is not sufficient

  l.prior <- 2*sqrt(beta/alpha/n0)*qt((1+level)/2,2*alpha)
  ifelse(l.prior <= len, 0, ceiling(-n0+qt((1+level)/2,2*alpha)^2*4*beta/alpha/len^2))
}

