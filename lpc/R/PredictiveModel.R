PredictiveModel <- function(dat, scores, thresh, type="regression"){
  ranked <- order(abs(scores), decreasing=TRUE)
  dat$x <- dat$x[ranked,]
  if(type=="survival"){
    p <- abs(cox.func(dat$x[1:thresh,], dat$y, dat$censoring.status, .05)$tt)
  } else if (type=="regression"){
    p <- abs(quantitative.func(dat$x[1:thresh,], dat$y, .05)$tt)
  }else if (type=="two class"){
    p <- abs(ttest.func(dat$x[1:thresh,], dat$y, .05)$tt)
  } else if (type=="multiclass"){
    p <- abs(multiclass.func(dat$x[1:thresh,], dat$y, .05)$tt)
  }
  pmeans <- NULL
  for(i in 1:length(p)) pmeans <- c(pmeans, mean(p[1:i]))
  return(pmeans)
}
