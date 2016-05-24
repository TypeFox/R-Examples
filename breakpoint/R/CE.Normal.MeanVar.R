CE.Normal.MeanVar <-
function(data, Nmax=10, eps=0.01, rho=0.05, M=200, h=5, a=0.8, b=0.8, distyp = 1, penalty = "BIC", parallel=FALSE){
  
  if(is.data.frame(data) == "FALSE"| is.null(dim(data)[2])) {
    print("Error in data : dataframe only")                                                           
  } else if(dim(data)[2] != 1) {
    print("Error in data : single column dataframe only")                               
  } else { 
    
  if(distyp == 1 & penalty == "BIC"){
      Melite <- M * rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.sim4beta.MeanVar.BIC", "betarand", "fun.alpha", "fun.beta", "llhood.MeanVarNormal", "loglik.MeanVarNormal", "BIC.MeanVarNormal"), envir=environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanVar.BIC(k, data, h, L0, L, M, Melite, eps, a)
            stopCluster(cl)   
            
          
          } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanVar.BIC(k, data, h, L0, L, M, Melite, eps, a)           
            
            
        } else sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta.MeanVar.BIC(k, data, h, L0, L, M, Melite, eps, a)
           
    BIC.summary <- sapply(sim, "[[", 2)
    opt.loci <- which(BIC.summary == min(BIC.summary))
      
    loci.BIC <- sim[[opt.loci]]$loci
    
    if(length(loci.BIC) >= 3) {
      return(list("No.BPs" = length(loci.BIC) - 2, "BP.Loc" = loci.BIC[2:(length(loci.BIC) - 1)], "BIC value" = sim[[opt.loci]]$BIC.Val, "ll" = sim[[opt.loci]]$LogLike))
    } else {
      return(paste("No Break-Points are Estimated")) 
      }
    
    } else if(distyp == 2 & penalty == "BIC"){
      
      Melite <- M*rho
      L <- length(data[, 1])
      L0 <- 1      
      k <- seq(0, Nmax, 1)
      
      if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.simnormal.MeanVar.BIC", "normrand", "llhood.MeanVarNormal", "loglik.MeanVarNormal", "BIC.MeanVarNormal"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b"), envir = environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanVar.BIC(k, data, h, L0, L, M, Melite, eps, a, b)
            stopCluster(cl)   
      
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanVar.BIC(k, data, h, L0, L, M, Melite, eps, a, b)
          
          
        } else { 
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal.MeanVar.BIC(k, data, h, L0, L, M, Melite, eps, a, b)
        }      
      
      BIC.summary <- sapply(sim, "[[", 2)
      opt.loci <- which(BIC.summary == min(BIC.summary))
      
      loci.BIC <- sim[[opt.loci]]$loci
      
        if(length(loci.BIC) >= 3) {
          return(list("No.BPs" = length(loci.BIC) - 2, "BP.Loc" = loci.BIC[2:(length(loci.BIC) - 1)], "BIC value" = sim[[opt.loci]]$BIC.Val, "ll" = sim[[opt.loci]]$LogLike))
        } else {
          return(paste("No Break-Points are Estimated")) 
        }          
    
      } else if (distyp == 1 & penalty == "AIC"){
        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- seq(0, Nmax, 1)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.sim4beta.MeanVar.AIC", "betarand", "fun.alpha", "fun.beta", "llhood.MeanVarNormal", "loglik.MeanVarNormal", "AIC.MeanVarNormal"), envir=environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanVar.AIC(k, data, h, L0, L, M, Melite, eps, a)
            stopCluster(cl)   
          
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
          registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.MeanVar.AIC(k, data, h, L0, L, M, Melite, eps, a)                           
          
          
        } else sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta.MeanVar.AIC(k, data, h, L0, L, M, Melite, eps, a)
        
        AIC.summary <- sapply(sim, "[[", 2)
        opt.loci <- which(AIC.summary == min(AIC.summary))

        loci.AIC <- sim[[opt.loci]]$loci
        
        if(length(loci.AIC) >= 3) {
          return(list("No.BPs" = length(loci.AIC) - 2, "BP.Loc" = loci.AIC[2:(length(loci.AIC) - 1)], "AIC value" = sim[[opt.loci]]$AIC.Val, "ll" = sim[[opt.loci]]$LogLike))
        } else {
          return(paste("No Break-Points are Estimated")) 
        }
        
      } else if (distyp == 2 & penalty == "AIC"){
        Melite <- M*rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- seq(0, Nmax, 1)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
            clusterExport(cl, c("ce.simnormal.MeanVar.AIC", "normrand", "llhood.MeanVarNormal", "loglik.MeanVarNormal", "AIC.MeanVarNormal"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b"), envir = environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanVar.AIC(k, data, h, L0, L, M, Melite, eps, a, b)
            stopCluster(cl)   
          
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
          registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.MeanVar.AIC(k, data, h, L0, L, M, Melite, eps, a, b)
          
          
        } else { 
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal.MeanVar.AIC(k, data, h, L0, L, M, Melite, eps, a, b)
        }      
        
        AIC.summary <- sapply(sim, "[[", 2)
        opt.loci <- which(AIC.summary == min(AIC.summary))
        
        loci.AIC <- sim[[opt.loci]]$loci
        
        if(length(loci.AIC) >= 3) {
          return(list("No.BPs" = length(loci.AIC) - 2, "BP.Loc" = loci.AIC[2:(length(loci.AIC) - 1)], "AIC value" = sim[[opt.loci]]$AIC.Val, "ll" = sim[[opt.loci]]$LogLike))
        } else {
          return(paste("No Break-Points are Estimated")) 
        }
  }
  }
}
