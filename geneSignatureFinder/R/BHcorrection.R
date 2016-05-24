BHcorrection <-
function(pvs, alpha = 0.05) {
    m <- length(pvs)
    hatK <- rev(which(m*pvs[order(pvs)]/1:m < alpha))[1]
    hatK <- ifelse(is.na(hatK), 1, hatK)
    ans <- pvs * m/hatK
    ans <- ifelse(ans >= 1, 1, ans)   
    attr(ans, "hatK") <- hatK
    attr(ans, "alpha") <- alpha
    return(ans)
  }
