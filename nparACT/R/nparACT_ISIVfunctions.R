nparACT_ISIVfunctions <- list(
  nparACT_ISIV = function(data_hrs, bin_hr){
    n <- nrow(data_hrs)
    p <- 1440/bin_hr  
    l <- 60/bin_hr   
    mean_all <- mean(data_hrs[1:n,]) 
    
    return_IS <- nparACT_ISIVfunctions$nparACT_IS(data_hrs, mean_all, bin_hr)
    
    IS <- return_IS[1]
    n <-  return_IS[3]
    p <-  return_IS[4]
    IV <- nparACT_ISIVfunctions$nparACT_IV(n, data_hrs, mean_all)
    
    result_ISIV <- c(IS, IV)
    return(result_ISIV)
  },
  
  nparACT_IS = function(data_hrs, mean_all, bin_hr){
    ## ---- IS numerator calculation
    result_ISnum <- matrix(NA, nrow = 24) 
    n <- nrow(data_hrs)
    p <- 1440/bin_hr  
    for (h in 1:24){ 
      s <- ceiling(n/p) 
      data_hrs3 <- data_hrs
      data_hrs3[s*p] <- NA 
      data_hrs3 <- matrix(data_hrs3) 
      hrlydat <- data_hrs3[c(seq(h,nrow(data_hrs3),24)),]
      hrlymean <- mean(hrlydat, na.rm = T) 
      x <- (hrlymean-mean_all)^2 
      result_ISnum[h,] <- x 
    }
    ISnum <- sum(result_ISnum)  
    ISnumerator <- n*ISnum  
        ## ---- IS denominator calculation
    result_ISdenom <- matrix(NA, nrow = n)  
    for (j in 1:n){
      y <- ((data_hrs[j,]-mean_all)^2)
      result_ISdenom[j,] <- y  
    }
    ISdenom <- sum(result_ISdenom)   
    ISdenominator <- p*ISdenom  
    ## -----------------------------
    IS <- round(ISnumerator/ISdenominator, digits = 2)
    return_IS <- c(IS, ISdenom, n, p)
    return(return_IS)
  },
  
  nparACT_IV = function(n, data_hrs, mean_all){
    result_IVnum <- matrix(NA, nrow = n) 
    for (k in 2:n){
      z <- ((data_hrs[k,]-data_hrs[(k-1),])^2)  
      result_IVnum[k,] <- z  
    }
    IVnum <- sum(result_IVnum, na.rm = T)  
    IVnumerator <- n*IVnum  
        ## ---- IV denominator calculation
    result_ISdenom <- matrix(NA, nrow = n)  
    for (j in 1:n){
      y <- ((data_hrs[j,]-mean_all)^2)
      result_ISdenom[j,] <- y  
    }
    ISdenom <- sum(result_ISdenom)   
    IVdenominator <- (n-1)*ISdenom ## ISdenom can be used!
    ## -----------------------------
    IV <- round(IVnumerator/IVdenominator, digits = 2)
    return(IV)
  }
)