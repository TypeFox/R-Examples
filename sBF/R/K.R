K <-
function(u, method = "gaussian") {
  switch(method,
    uniform=g<-0,
    epanechnikov=g<-1,
    biweight=g<-2,
    triweight=g<-3,
    return(dnorm(u))
  )
  return(1/beta(1/2, g+1)*(1-u^2)^g*(abs(u)<=1))
}

