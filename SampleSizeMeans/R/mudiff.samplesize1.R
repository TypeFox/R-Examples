`mudiff.samplesize1` <-
function(alpha1,beta1,alpha2,beta2,n01,n02,equal,n1)
{
  n2 <- ifelse(equal, n1, sqrt(beta2/beta1*(alpha1-1)/(alpha2-1))*(n1+n01)-n02)

  max(0,ceiling(n2))
}

