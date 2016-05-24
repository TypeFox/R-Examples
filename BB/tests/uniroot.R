
# This test will fail without try() wrapper in multiStart()
require(BB)

test <- function(x, bb0=-3, bb1=5, c0=2, r0=0) { 
  ((exp(c0-r0)*(bb0+x)*(bb1-x))/((bb0+x+1)*(bb1-x-1))-1)
  } 

# uniroot(test,c(-100,100))$root   this will fail

p0 <- matrix(seq(0, 6, length=10), ncol=1)
ans   <- multiStart(par=p0, fn=test, method=2, control=list(M=20, NM=FALSE)) 
ans$par[ans$conv]

