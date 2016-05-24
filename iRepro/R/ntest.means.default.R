ntest.means.default <-
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
  
  t <- .chi2.int.test(rowMeans(avg.lower), rowMeans(avg.upper), bins, est$mu, sqrt(est$sigma2.b + 0.5*est$sigma2.w))  
  test.res <- list(statistic = t$statistic,
                  parameter = t$parameter,
                  p.value = t$p.value,
                  data = "means",
                  mu = est$mu,
                  var = est$sigma2.b + 0.5*est$sigma2.w,
                  bins = bins)
  class(test.res) <- "ntestMeans"  
  test.res
}

