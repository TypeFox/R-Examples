truncSum <- function(stim, y, hrf, isinterp=F, iscor=F) {
  yn <- length(y)
  ns <- length(stim)

  agg <- rep(0, ns)
  if(!isinterp) stim <- round(stim)

  nh1 <- length(hrf)-1

  for (i in 1:ns) {
    if (isinterp) {
      tmp <- stim[i]:(stim[i]+nh1)
      tmp <- tmp[tmp<=yn]
      ys <- spline(1:yn,y,xout=tmp)$y
    } else {
      ys <- y[stim[i]:min(yn, (stim[i]+nh1))]
    }
    ysn <- length(ys)
    hr <- hrf[1:ysn]
    if (iscor) {
      tmp <- sqrt(sum(hr^2)*sum(ys^2))
      if(tmp==0) {
        agg[i] <- 0
      } else {
        agg[i] <- sum(hr*ys)/tmp
      }
    } else {
      if (sum(hr)==0) {
        hr <- 0*hr
      } else {
        hr <- hr/sum(hr)
      }
      agg[i] <- sum(hr*ys)
    }
  }
  return(agg)
}
