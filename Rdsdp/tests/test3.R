# Test a simple example
library(Rdsdp)

    # Sedumi format example
    K=c()
    K$s=c(3)
    K$l=0
    OPTIONS=NULL
    OPTIONS$gaptol=1e-10

  A1 = c(2,-0.5,-0.6, -0.5,2,0.4, -0.6,0.4,3)
  C = -A1
  # A1 = c(0,1,0, 1,0,0, 0,0,0)
  A2 = c(0,0,1, 0,0,0, 1,0,0)
  A3 = c(0,0,0, 0,0,1, 0,1,0)
  A4 = -c(1,0,0, 0,1,0, 0,0,1)
  A=rbind(A2,A3,A4)
  b <- -c(0,0,1)

  ret = dsdp(A,b,C,K,OPTIONS)
  y = ret$y

  stopifnot(all.equal(y,c(0.6,-0.4,3),tolerance=1e-05))



