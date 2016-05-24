plotbe <- function(m, mpl, return.be = FALSE, ...){
  baseline <- basehaz(m, centered = FALSE)
  plot(baseline$time, baseline$hazard, type = "n", xlab = "Time", bty = "n",
       ylab = "Cumulative Baseline Hazard", ...)
  for(x.q in 1:length(levels(baseline$strata))){
    time <- baseline$time[baseline$strata == levels(baseline$strata)[x.q]]
    cbhr <- baseline$hazard[baseline$strata == levels(baseline$strata)[x.q]]
    lines(time, cbhr, type = "s", lty = x.q)
  }
  ## same legend labels as in nae:
  k <- l <- q <- NULL # rep(NA, length(mpl)*length(mpl))
  for(k.i in 1:length(mpl)){
    all.to <- as.vector(mpl[[k.i]]$all.to)
    if(!is.null(all.to)){
      for(l.i in 1:length(all.to)){
        k <- c(k, k.i)
        l <- c(l, all.to[l.i])
        q <- c(q, paste(k.i, all.to[l.i], sep = ""))
      }
    }
  }
  legend.text <- unlist(strsplit(levels(baseline$strata), 
                                 split = "="))[seq(2, 2*length(levels(baseline$strata)), by = 2)]
  for(q.i in 1:length(legend.text)){
    hi <- which(legend.text[q.i] == q)
    legend.text[q.i] <- paste(k[hi], l[hi], sep = " ")
  }
  legend("topleft", lty = 1:length(levels(baseline$strata)), 
         legend.text, bty = "n")
  if(return.be){
    return(baseline)
  }
}