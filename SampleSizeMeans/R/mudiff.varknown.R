`mudiff.varknown` <-
function(len,lambda1,n01,lambda2,n02,level=0.95,equal=TRUE)
{
  l.prior <- 2*qnorm((1+level)/2)*sqrt(1/n01/lambda1+1/n02/lambda2)
  
  if (l.prior <= len)
  {
    0 # prior knowledge is sufficient
  }
  else
  {
    lambda01 <- n01*lambda1
    lambda02 <- n02*lambda2

    d <- len*len/qnorm((level+1)/2)^2/4
    a <- lambda1*lambda2*d
    b <- (lambda1*lambda02+lambda2*lambda01)*d-lambda1-lambda2
    cc <- lambda01*lambda02*d-lambda01-lambda02
    common.n <- ceiling((-b+sqrt(b*b-4*a*cc))/2/a)

    if (equal)
    { n1 <- common.n
      n2 <- common.n
    }
    else
    {
      n1 <- (1+sqrt(lambda1/lambda2))/d/lambda1-n01
      n1 <- max(0,n1)
      n2 <- sqrt(lambda1/lambda2)*(n01+n1)-n02
      n2 <- max(0,n2)

      # The optimal solution is on one of the three upper corners of square
      # around (n1,n2) [if n1 and n2 are not integers]

      if (n1!=floor(n1) || n2!=floor(n2))
      {
        vn1 <- c(floor(n1),ceiling(n1),ceiling(n1))
        vn2 <- c(ceiling(n2),floor(n2),ceiling(n2))

        v <- 1/(lambda01+lambda1*vn1)+1/(lambda02+lambda2*vn2)

        n1 <- vn1[v<=d][1]
        n2 <- vn2[v<=d][1]
      }

      if (n2==0) n1 <- ceiling((lambda01+lambda02-d*lambda01*lambda02)/lambda1/(d*lambda02-1))
      if (n1==0) n2 <- ceiling((lambda01+lambda02-d*lambda01*lambda02)/lambda2/(d*lambda01-1))
    }

    c(n1,n2)
  }
}

