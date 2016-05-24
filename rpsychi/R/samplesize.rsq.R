samplesize.rsq <- function(delta, n.ind, power = .80, sig.level=.05){
    fsq <- delta/(1 - delta)
    temp <- function(n){
      v <- n - n.ind - 1
      
      lamda <- fsq * (n.ind + v + 1)
      return(pf(qf(p=sig.level, df1=n.ind, df2=v, lower.tail=FALSE),df1=n.ind, df2=v, ncp=lamda, lower.tail=FALSE) - power)
      }
    output <- ceiling(uniroot(temp,c(n.ind+2, 10e10))[[1]])
    return(output)
}

