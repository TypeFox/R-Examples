library(lgcp)

tfuncs <- function(n,f,s=10){
  x = 1:s; y = x + 2

  for(i in 1:n){
    g = list(x=x,y=y,z=matrix(runif(s*s),s,s))
    f(g,i)
  }

  return(f)
  
}

tstt <- function(n,s,seed=1){
  set.seed(seed)

  testList = list(sumsumsq(),keepall(),exceed(c(0.1,0.2,0.3)),minmax())

  for(t in testList){
    cat("testing...\n")
    print(t)
    t1 = tfuncs(n,t,s=s)
    print(getResults(t1))
  }
}
