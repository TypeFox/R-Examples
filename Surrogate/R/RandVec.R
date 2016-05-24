
RandVec <- function(a=0, b=1, s=1, n=9, m=1,  Seed=sample(1:1000, size=1)) {    
  
  # This function is an R adaptation of a matlab program written by Roger Stafford 
  # For details on the original Matlab algorithm, see: 
  # http://www.mathworks.com/matlabcentral/fileexchange/9700-random-vectors-with-fixed-sum/content/randfixedsum.m
  
  set.seed(Seed)
  
  if (b-a==0){x <- matrix(rep(0, times=(n*m)), ncol=m)}  #added
  
  if (b-a!=0){
    s <- (s-n*a)/(b-a)
    
    k <- max(min(floor(s), n-1), 0)
    s <- max(min(s, k+1), k)
    s1 <- s - seq(from=k, to=(k-n+1), by=-1)
    s2 <- seq(from=(k+n), to=(k+1), by=-1) - s
    w <- matrix(data=(rep(0, times=(n*(n+1)))), ncol=(n+1), nrow=n)
    realmax <- 1e+300 #largest integer; modified
    w[1, 2] <- realmax
    t <- matrix(data=(rep(0, times=((n-1)*(n)))), ncol=(n), nrow=(n-1))
    tiny <- 2^(-1074)
    
    for (i in 2:n) {
      tmp1 <- (w[i-1,])[2:(i+1)]*s1[1:i]/i
      tmp2 <- (w[i-1,])[1:i]*s2[(n-i+1):n]/i
      tmp <- as.vector(tmp1 + tmp2)
      w[i,][2:(i+1)] <- tmp 
      tmp3 <- w[i,][2:(i+1)] + tiny
      tmp4 <- as.numeric(s2[(n-i+1):n] > s1[1:i])
      t[i-1,][1:i] <- (tmp2/tmp3)*tmp4 + (1-tmp1/tmp3)*rev(tmp4)
    }
    
    v <- n^(3/2)*((w[n, k+2])/realmax)*(b-a)^(n-1)
    x <- matrix(data=(rep(0, times=(n*m))), ncol=m, nrow=n)
    if (m<1) {stop("m should be >= 1")}
    rt <- matrix(runif(m*(n-1)), ncol=m)
    rs <- matrix(runif(m*(n-1)), ncol=m)
    s <- rep(s, times=m)
    j <- rep(k+1, times=m)
    sm <- matrix(data=(rep(0, times=(1*m))), ncol=m, nrow=1)
    pr <- matrix(data=(rep(1, times=(1*m))), ncol=m, nrow=1)
    
    for (i in (n-1):1) {
      e <- as.numeric(rt[n-i,] <=  rep(t[i,j[1]], times=length(j)))
      sx <- rs[n-i,]^(1/i)
      sm <- sm + (1-sx)*pr*s/(i+1)
      pr = sx*pr 
      x[n-i,] <- sm + pr*e
      s <- s-e
      j <- j-e
    }
    
    x[n,] <- sm + pr*s
    
    # permute and rescale
    set.seed(Seed)
    rp <- matrix(runif(m*n), ncol=m)
    offset <- rep(seq_len(m), rep(n, m))
    rp <- matrix(sort(rp + offset), nrow=nrow(rp)) - offset 
    
    p <- matrix(sample(1:n, n, replace=F), ncol=1)
    if (m>1){
      for (i in 1:(m-1)){
        p <- cbind(p, matrix(sample(1:n, n,replace=F), ncol=1))
      }
    }
    
    repmat <- seq(from=0, to=(n*(m-1)), by=n)
    for (i in 1:(n-1)){
      repmat <- rbind(seq(from=0, to=(n*(m-1)), by=n), repmat)
    }
    repmat <- matrix(repmat, ncol=m)
    
    x <- (b-a)*matrix(x[order(p+repmat)], ncol=m)+a
  }
  
  
  fit <- list(RandVecOutput=x)   
  
  class(fit) <- "RandVec"
  fit
}
