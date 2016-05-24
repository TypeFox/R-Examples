efaTable <- 
  function(fa.res, data, mean.sd = TRUE) {
    tab <- round(data.frame(fa.res$loadings[,]), 2)
    tab <- data.frame(tab, h2 = round(fa.res$communality, 2))
    
    if(mean.sd == TRUE) {
      data.m <- plyr::ldply(data, mean, na.rm=TRUE)
      data.sd <- plyr::ldply(data, stats::sd, na.rm=TRUE)
      
      tab <- data.frame(M = round(data.m[, 2], 2), SD = round(data.sd[, 2], 2), tab)
    }
    class(tab) <- c("data.frame", paste0("efaTable", ifelse(mean.sd, ".m.sd", "")), paste0(ifelse(mean.sd, ncol(tab) - 3, ncol(tab) - 1), "factors"))
    
    tab
    
  }
