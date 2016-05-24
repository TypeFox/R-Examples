# estimator based on
# Antonius M. J. VanDongen (1996) "A New Algorithm for Idealizing Single Ion Channel Data Containing Multiple Unknown Conductance Levels", Biophysical Journal 70, 1303-1315.
"transit" <-
function(y, x = 1:length(y), x0 = 2 * x[1] - x[2], sigma.amp = NA, sigma.slope = NA, amp.thresh = 3, slope.thresh = 2, rel.amp.n = 3, rel.amp.thresh = 4, family = c("gauss", "gaussKern"), param = NULL, refit = FALSE)
{
  family <- match.arg(family)
  
  # compute candidates (for output later)
  cand <- stepcand(y = y, x = x, x0 = x0, family = family, param = param)

  # length of kernel (if any)
  kl <- if(family == "gaussKern") length(param$kern) else 1
  
  # central difference
  cdif <- diff(y, lag = 2) / 2
  
  # if no sigma has been specified use sdrobnorm
  if(is.na(sigma.amp)) sigma.amp <- sdrobnorm(y,, kl)
  if(is.na(sigma.slope)) sigma.slope <- sdrobnorm(cdif,, kl + 2)
  
  # pad cdif with 0s to length of y
  cdif <- c(0, cdif, 0)
  
  # critical cdif (Fig. 3B)
  crit <- abs(cdif) > slope.thresh * sigma.slope
  
  # no forward difference used as VanDongen does not explain precisely how...
  
  # determine left and right indices of levels
  if(any(crit)) {
    li <- which(!crit & c(TRUE, crit[-length(y)]))
    ri <- which(!crit & c(crit[-1], TRUE))
  } else {
    li <- 1
    ri <- length(y)
  }
  
  # value is obtained by averaging over levels
  value <- sapply(1:length(li), function(j) mean(y[li[j]:ri[j]]))
  len <- ri - li + 1
  
  # remove spurious transitions
  while(length(ri) > 1) {
    # the minimal difference in value divided by sigma.amp and amp.tresh (or rel.amp.tresh if one level is shorter than rel.amp.n)
    dvalue <- abs(diff(value)) / sigma.amp / ifelse(len[-length(len)] <= rel.amp.n | len[-1] <= rel.amp.n, rel.amp.thresh, amp.thresh)
    ind <- which.min(dvalue)[1]
    if(dvalue[ind] < 1) {
      value[ind + 1] <- ( value[ind] * len[ind] + value[ind + 1] * len[ind + 1] ) / ( len[ind] + len[ind + 1] )
      value <- value[-ind]
      len[ind + 1] <- len[ind] + len[ind + 1] # while not written in VanDongen, appears consequent since he does not compute the new value over the entire new level either
      len <- len[-ind]
      li[ind + 1] <- li[ind]
      li <- li[-ind]
      ri <- ri[-ind]
    } else break
  }
  
  # add half a transition's length to each neighbouring level
  ri[-length(ri)] <- ri[-length(ri)] + ( li[-1] - ri[-length(ri)] ) %/% 2

  # ouput
  ret <- stepfit(cost = NA, family = attr(cand, "family"), value = value, param = attr(cand, "param"), 
    leftEnd = cand$leftEnd[c(0, ri[-length(ri)]) + 1], rightEnd = cand$rightEnd[ri], x0 = attr(cand, "x0"),
    leftIndex = cand$leftIndex[c(0, ri[-length(ri)]) + 1], rightIndex = cand$rightIndex[ri])
  if(attr(cand, "family") != "gaussvar") {
    ret$cumSum <- cand$cumSum[ri]
  }
  ret$cumSumWe <- cand$cumSumWe[ri]
  if(attr(cand, "family") %in% c("gauss", "gaussvar")) {
    ret$cumSumSq <- cand$cumSumSq[ri]
  }
  if(attr(cand, "family") == "gaussKern") {
    ret$cumSumSq <- cand$cumSumSq[ri]
    ret$lXy <- cand$lXy[ri]
    ret$lcXy <- cand$lcXy[ri]
    ret$rcXy <- cand$rcXy[ri]
    ret$rXy <- cand$rXy[ri]
    if(!identical(refit, FALSE)) ret <- ret[refit = refit]
  }
  ret
}
