
lpc.curve.length <- function(lpcsl,or.pi,branch=0,total.subdivisions=1e4,min.subdivisions=100) {
    gradnorm <- function(tt,sps) {
      sum <- numeric(length(tt))
      for (j in 1:length(sps))
        sum <- sum + sps[[j]](tt,deriv=1)^2
      sqrt(sum)
    }
    if (length(branch)==1)
      branch <- rep(branch,length(or.pi))
    if (length(or.pi)!=length(branch))
      stop("Length mismatch between or.pi and branch")
    result <- numeric(length(or.pi))
    branches <- unique(branch)
    for (cur.branch in branches) {
      cur.or.pi <- or.pi[branch==cur.branch]
      length.between.points <- numeric(length(cur.or.pi))
      cur.or.pi.order <- order(cur.or.pi)
      cur.or.pi <- cur.or.pi[cur.or.pi.order]
      from <- c(0,cur.or.pi)
      for (i in 1:length(cur.or.pi.order)) {        
        subdivisions <- max(min.subdivisions,ceiling(total.subdivisions*(cur.or.pi[i]-from[i])/(lpcsl[[cur.branch+1]]$range[2]-lpcsl[[cur.branch+1]]$range[1])))
        if(from[i]==cur.or.pi[i]){
                length.between.points[i] <-0
        } else {    
            length.between.points[i] <- integrate(gradnorm,lower=from[i],upper=cur.or.pi[i],subdivisions=subdivisions,sps=lpcsl[[cur.branch+1]]$splinefun)$value
        }
      }
      total.length <- cumsum(length.between.points)
      result[branch==cur.branch][cur.or.pi.order] <- total.length
    }
    result
  } 

