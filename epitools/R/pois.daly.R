"pois.daly" <-
  function(x, pt = 1, conf.level = 0.95){
    xc <- cbind(x,conf.level,pt)
    pt2 <- xc[,3]
    results <- matrix(NA,nrow(xc),6)
    cipois <-
      function(x, conf.level = 0.95){
        if(x!=0){
          LL <- qgamma((1 - conf.level)/2, x)
          UL <- qgamma((1 + conf.level)/2, x + 1)
        } else {
          if(x==0){
            LL <- 0
            UL <- -log(1 - conf.level)
          }
        }
        data.frame(x = x, lower = LL, upper = UL)
      }
    for(i in 1:nrow(xc)){
      alp <- 1-xc[i,2]
      daly <- cipois(x = xc[i, 1], conf.level = xc[i, 2])
      LCL <- daly$lower/pt2[i]
      UCL <- daly$upper/pt2[i]
      results[i,] <- c(xc[i,1],pt2[i],xc[i,1]/pt2[i],LCL,UCL,xc[i,2]) 
    }
    coln <- c("x","pt","rate","lower","upper","conf.level")
    colnames(results) <- coln
    data.frame(results)
  }
