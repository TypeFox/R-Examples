CE.Normal.Mean <-
function(data = data, Nmax=10, eps=0.01, rho=0.05, M=200, h=5, a=0.8, b=0.8, distyp = 1, penalty = "mBIC", parallel=FALSE){
  
  if(is.data.frame(data) == "FALSE"| is.null(dim(data)[2])) {
    print("Error in data : dataframe only")                                                           
  } else if(dim(data)[2] != 1) {
    print("Error in data : single column dataframe only")                               
  } else { 
    
    if(distyp == 1 & penalty == "mBIC"){
      Melite <- M * rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.sim4beta", "betarand", "fun.alpha", "fun.beta", "mBIC"), envir=environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta(k, data, h, L0, L, M, Melite, eps, a)
            stopCluster(cl)   
            
          
          } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta(k, data, h, L0, L, M, Melite, eps, a)        
            
        } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta(k, data, h, L0, L, M, Melite, eps, a)
        } 
    mBic.summary <- sapply(sim, "[[", 2)
    opt.loci <- which(mBic.summary == max(mBic.summary))
    
    loci.mBIC <- sim[[opt.loci]]$loci

    logLL <- llhood.MeanNormal(loci.mBIC, data, v=var(data[ ,1]) , h)

    if(length(loci.mBIC) >= 3) {
      return(list("No.BPs" = length(loci.mBIC) - 2, "BP.Loc" = loci.mBIC[2:(length(loci.mBIC) - 1)], "mBIC" = sim[[opt.loci]]$mBIC , "ll" = logLL))
    } else {
      return(paste("No Break-Points are Estimated")) 
      }
    
    } else if(distyp == 2 & penalty == "mBIC"){
      
      Melite <- M*rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
      if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.simnormal", "normrand", "mBIC"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b"), envir = environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal(k, data, h, L0, L, M, Melite, eps, a, b)
            stopCluster(cl)   
          
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal(k, data, h, L0, L, M, Melite, eps, a, b)
          
          
        } else { 
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal(k, data, h, L0, L, M, Melite, eps, a, b)
        }      
      
      mBic.summary <- sapply(sim, "[[", 2)
      opt.loci <- which(mBic.summary == max(mBic.summary))
      
      loci.mBIC <- sim[[opt.loci]]$loci
      logLL <- llhood.MeanNormal(loci.mBIC, data, v=var(data[ ,1]) , h)
      
        if(length(loci.mBIC) >= 3) {
          return(list("No.BPs" = length(loci.mBIC) - 2, "BP.Loc" = loci.mBIC[2:(length(loci.mBIC) - 1)], "mBIC" = sim[[opt.loci]]$mBIC, "ll" = logLL))
        } else {
          return(paste("No Break-Points are Estimated")) 
        }          
    
      } else if(distyp == 1 & penalty == "BIC") {
        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- seq(0, Nmax, 1)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.sim4beta.MeanBIC", "betarand", "fun.alpha", "fun.beta", "llhood.MeanNormal", "loglik.MeanNormal", "BIC.MeanNormal"), envir=environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanBIC(k, data, h, L0, L, M, Melite, eps, a)
            stopCluster(cl)   
          
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanBIC(k, data, h, L0, L, M, Melite, eps, a)
          
        } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta.MeanBIC(k, data, h, L0, L, M, Melite, eps, a)
        }
        
        Bic.summary <- sapply(sim, "[[", 2)
        opt.loci <- which(Bic.summary == min(Bic.summary))

        loci.BIC <- sim[[opt.loci]]$loci

        if(length(loci.BIC) >= 3) {
          return(list("No.BPs" = length(loci.BIC) - 2, "BP.Loc" = loci.BIC[2:(length(loci.BIC) - 1)], "BIC" = sim[[opt.loci]]$BIC.Val, "ll" = sim[[opt.loci]]$LogLike ))
        } else {
          return(paste("No Break-Points are Estimated")) 
        }
        
      } else if(distyp == 2 & penalty == "BIC"){
        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- seq(0, Nmax, 1)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.simnormal.MeanBIC", "normrand", "llhood.MeanNormal", "loglik.MeanNormal", "BIC.MeanNormal"), envir=environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanBIC(k, data, h, L0, L, M, Melite, eps, a, b)
            stopCluster(cl)   
          
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanBIC(k, data, h, L0, L, M, Melite, eps, a, b)                           
          
          
        } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal.MeanBIC(k, data, h, L0, L, M, Melite, eps, a, b)
        }
        
        Bic.summary <- sapply(sim, "[[", 2)
        opt.loci <- which(Bic.summary == min(Bic.summary))
        
        loci.BIC <- sim[[opt.loci]]$loci
        
        if(length(loci.BIC) >= 3) {
          return(list("No.BPs" = length(loci.BIC) - 2, "BP.Loc" = loci.BIC[2:(length(loci.BIC) - 1)], "BIC" = sim[[opt.loci]]$BIC.Val, "ll" = sim[[opt.loci]]$LogLike ))
        } else {
          return(paste("No Break-Points are Estimated")) 
        }
        
    } else  if(distyp == 1 & penalty == "AIC"){
      Melite <- M * rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
      if(parallel == TRUE & .Platform$OS.type == "windows"){
          cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.sim4beta.MeanAIC", "betarand", "fun.alpha", "fun.beta", "llhood.MeanNormal", "loglik.MeanNormal", "AIC.MeanNormal"), envir=environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanAIC(k, data, h, L0, L, M, Melite, eps, a)
          stopCluster(cl)   

      } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanAIC(k, data, h, L0, L, M, Melite, eps, a)          
      } else {
        sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta.MeanAIC(k, data, h, L0, L, M, Melite, eps, a)
      }
      Aic.summary <- sapply(sim, "[[", 2)
      opt.loci <- which(Aic.summary == min(Aic.summary))
      loci.AIC <- sim[[opt.loci]]$loci
      
      if(length(loci.AIC) >= 3) {
        return(list("No.BPs" = length(loci.AIC) - 2, "BP.Loc" = loci.AIC[2:(length(loci.AIC) - 1)], "AIC" = sim[[opt.loci]]$AIC.Val, "ll" = sim[[opt.loci]]$LogLike ))
      } else {
        return(paste("No Break-Points are Estimated")) 
      }
      
    } else  if(distyp == 2 & penalty == "AIC"){
      Melite <- M * rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
      if(parallel == TRUE & .Platform$OS.type == "windows"){
          cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.simnormal.MeanAIC", "normrand", "llhood.MeanNormal", "loglik.MeanNormal", "AIC.MeanNormal"), envir=environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanAIC(k, data, h, L0, L, M, Melite, eps, a, b)
          stopCluster(cl)   
        
        
      } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanAIC(k, data, h, L0, L, M, Melite, eps, a, b)                         
        
      } else sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal.MeanAIC(k, data, h, L0, L, M, Melite, eps, a, b)
      
      Aic.summary <- sapply(sim, "[[", 2)
      opt.loci <- which(Aic.summary == min(Aic.summary))
      loci.AIC <- sim[[opt.loci]]$loci
      
      if(length(loci.AIC) >= 3) {
        return(list("No.BPs" = length(loci.AIC) - 2, "BP.Loc" = loci.AIC[2:(length(loci.AIC) - 1)], "AIC" = sim[[opt.loci]]$AIC.Val, "ll" = sim[[opt.loci]]$LogLike))
      } else {
        return(paste("No Break-Points are Estimated")) 
      }
    }
  }
}
