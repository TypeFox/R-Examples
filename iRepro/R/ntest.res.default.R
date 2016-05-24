ntest.res.default <-
function(r1, r2, predefined.classes=FALSE, classes, c.limits, optim.method=1, bins=10){
  
  bins <- round(bins)
  if(bins<0)
    stop("Number of bins must be a positive integer")
  
  est <- intervalICC(r1, r2, predefined.classes, classes, c.limits, optim.method)  
  ratings <- cbind(r1, r2)
  if(any(is.na(ratings))){    
    ratings <- na.omit(ratings)
  }
  
  if(predefined.classes){
    avg.lower <- ratings
    avg.upper <- ratings
    for(i in 1:length(classes)){
      avg.lower[avg.lower==classes[i]] <- c.limits[i,1]
      avg.upper[avg.upper==classes[i]] <- c.limits[i,2]
    }
  }else{
    avg.lower <- ratings[,c(1,3)]
    avg.upper <- ratings[,c(2,4)]
  }
  
  
  
  avg.l <- rowMeans(avg.lower)
  avg.u <- rowMeans(avg.upper)
  t1 <- .chi2.int.test(avg.lower[,1]-avg.u, avg.upper[,1]-avg.l, bins, 0, sqrt(0.5*est$sigma2.w))
  t2 <- .chi2.int.test(avg.lower[,2]-avg.u, avg.upper[,2]-avg.l, bins, 0, sqrt(0.5*est$sigma2.w))

  test.res <- list(statistic.res1 = t1$statistic,
                  p.value.res1 = t1$p.value,
                  statistic.res2 = t2$statistic,
                  p.value.res2 = t2$p.value,
                  parameter = t1$parameter,
                  data = "residuals",
                  mu = 0,
                  var = 0.5*est$sigma2.w,
                  bins = bins) 
  class(test.res) <- "ntestRes"
  test.res
}

