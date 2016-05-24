jointcrit=function (p1, p2, alpha = 0.05) {
  p <- gconv(p1, p2)
  rcp <- cumsum(rev(p))
  rcp1 <- cumsum(rev(p1))
  if((rcp[1]>alpha & rcp1[1]>alpha)){
    out=matrix(c(Inf,Inf,0,0,0,0),nrow=1)
    colnames(out) <- c("t1", "t1+t2", "Pr(T1>=t1)", "Pr(T1+T2>=t1+t2)", 
                       "Pr(T1>=t1 or T1+T2>=t1+t2)", "|Difference|")
  }
  if((rcp[1])<=alpha & rcp1[1]>alpha){
    no.rcp.leq.alpha=length(rcp[rcp<=alpha])
    kmax <- length(p) - 1
    kmin <- kmax + 1 - no.rcp.leq.alpha    
    out=matrix(c(Inf,kmin,0,rcp[no.rcp.leq.alpha],rcp[no.rcp.leq.alpha],rcp[no.rcp.leq.alpha]),nrow=1)
    colnames(out) <- c("t1", "t1+t2", "Pr(T1>=t1)", "Pr(T1+T2>=t1+t2)", 
                       "Pr(T1>=t1 or T1+T2>=t1+t2)", "|Difference|")
  }
  if((rcp[1])<=alpha & rcp1[1]<=alpha){
    kmax <- length(p) - 1
    kmin <- kmax + 1 - length(rcp[rcp <= alpha])
    jnt <- outer(p1, p2, "*")
    n1 <- length(p1) - 1
    n2 <- length(p2) - 1
    idx <- outer(0:n1, 0:n2, "+")
    idxj <- outer(0:n1, rep(1, n2 + 1), "*")
    jmin <- n1 + 1 - length(rcp1[rcp1 <= alpha])
    cp <- rev(rcp)
    cp1 <- rev(rcp1)
    out <- matrix(NA, (kmax + 1 - kmin) * (n1 + 1 - jmin), 
                  6)
    colnames(out) <- c("t1", "t1+t2", "Pr(T1>=t1)", "Pr(T1+T2>=t1+t2)", 
                       "Pr(T1>=t1 or T1+T2>=t1+t2)", "|Difference|")
    i <- 0
    for (k in kmin:kmax) {
      for (j in jmin:n1) {
        i <- i + 1
        out[i, 1] <- j
        out[i, 2] <- k
        out[i, 3] <- cp1[j + 1]
        out[i, 4] <- cp[k + 1]
        out[i, 5] <- sum(jnt[(idx >= k) | (idxj >= j)])
        out[i, 6] <- abs(out[i, 3] - out[i, 4])
      }
    }
    out <- out[out[, 5] <= alpha, ]
    out <- out[order(out[, 1] + out[, 2], out[, 1]), ]
  }
  out
}

