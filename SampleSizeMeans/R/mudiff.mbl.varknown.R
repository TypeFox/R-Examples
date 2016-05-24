`mudiff.mbl.varknown` <-
function(len,lambda1,lambda2,level=0.95,equal=TRUE)
{
  common.n <- ceiling(4*qnorm((1+level)/2)^2*(1/lambda1+1/lambda2)/len^2)

  if (equal) n <- c(common.n,common.n)
  else
  {
    n1 <- 4/len^2*qnorm((1+level)/2)^2*(1/lambda1+1/sqrt(lambda1*lambda2))
    n2 <- sqrt(lambda1/lambda2)*n1

    # The optimal choice is one of the three(up-left, up-right or down-right)
    # corners of the square surrounding (n1,n2)
    # [if n1 and n2 are not integers]

    if (n1!=floor(n1) || n2!=floor(n2))
    {
      n <- matrix(c(floor(n1),ceiling(n1),ceiling(n1),
                 ceiling(n2),floor(n2),ceiling(n2)),ncol=2)

      n <- n[2*qnorm((1+level)/2)*sqrt(1/lambda1/n[,1]+1/lambda2/n[,2])<=len]

      # If n is not a vector, it is a matrix and its first line
      # contains optimal sample sizes

      if (!is.vector(n)) n <- n[1,]
    }
  }

  n
}

