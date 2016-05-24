CE.ZINB.Init <-
function(data, init.locs, eps=0.01, rho=0.05, M=200, h=5, a=0.8, b=0.8, distyp = 1, penalty = "BIC", var.init = 100000, parallel=FALSE){
    
  if(is.data.frame(data) == "FALSE"| is.null(dim(data)[2])) {
    print("Error in data: dataframe only")                                                           
  } else if(dim(data)[2] != 1) {
    print("Error in data: single column dataframe only")                               
  } else if(missing(init.locs)){
    print("Error: Initial locations are not provided!!!")
  
    } else { 
      
      if(distyp == 1 & penalty == "BIC"){
        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        r <- suppressWarnings(try(fitdistr(data[, 1], "negative binomial")[[1]][[1]], silent = T)) 

        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type = "SOCK") 
            clusterExport(cl, c("ce.4betaZINB.BIC.Init", "betarand", "fun.alpha", "fun.beta", "BICzinb", "llhoodzinb", "loglikzinb", "betaIntEst"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "init.locs", "r", "var.init"), envir = environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.4betaZINB.BIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, r, var.init)
            stopCluster(cl)   
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.4betaZINB.BIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, r, var.init)
                                     
          
        } else sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.4betaZINB.BIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, r, var.init)
        
        cpt.loci <- sim[[1]]$loci
        
        return(list("No.BPs" = length(cpt.loci) - 2, "BP.Loc" = cpt.loci[2:(length(cpt.loci) - 1)], "BIC" = sim[[1]]$BIC, "ll" = sim[[1]]$LogLike))
        
      } else if(distyp == 2 & penalty == "BIC"){
        
        Melite <- M*rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        r <- suppressWarnings(try(fitdistr(data[, 1], "negative binomial")[[1]][[1]], silent = T)) 
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type = "SOCK") 
            clusterExport(cl, c("ce.simNormalZINB.BIC.Init", "normrand", "BICzinb", "llhoodzinb", "loglikzinb"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b",  "init.locs", "r", "var.init"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simNormalZINB.BIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, r, var.init)
            stopCluster(cl)

          } else if (parallel == TRUE & .Platform$OS.type == "unix"){
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simNormalZINB.BIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, r, var.init)
            
            
          } else {
            sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simNormalZINB.BIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, r, var.init)
          }
        
        cpt.loci <- sim[[1]]$loci
        
        return(list("No.BPs" = length(cpt.loci) - 2, "BP.Loc" = cpt.loci[2:(length(cpt.loci) - 1)], "BIC" = sim[[1]]$BIC, "ll" = sim[[1]]$LogLike))
        
      } else if(distyp == 1 & penalty == "AIC"){
      
        Melite <- M*rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        r <- suppressWarnings(try(fitdistr(data[, 1], "negative binomial")[[1]][[1]], silent = T)) 
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type = "SOCK") 
            clusterExport(cl, c("ce.4betaZINB.AIC.Init", "betarand", "fun.alpha", "fun.beta", "AICzinb", "llhoodzinb", "loglikzinb", "betaIntEst"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "init.locs", "r", "var.init"), envir = environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.4betaZINB.AIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, r, var.init)
            stopCluster(cl)   

        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
           registerDoParallel(parallel::detectCores()) 
           sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.4betaZINB.AIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, r, var.init)
          
        } else sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.4betaZINB.AIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, r, var.init)
        
        cpt.loci <- sim[[1]]$loci

        return(list("No.BPs" = length(cpt.loci) - 2, "BP.Loc" = cpt.loci[2:(length(cpt.loci) - 1)], "AIC" = sim[[1]]$AIC, "ll" = sim[[1]]$LogLike))
        
      } else if(distyp == 2 & penalty == "AIC"){
       
        Melite <- M*rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        r <- suppressWarnings(try(fitdistr(data[, 1], "negative binomial")[[1]][[1]], silent = T)) 
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type = "SOCK") 
            clusterExport(cl, c("ce.simNormalZINB.AIC.Init", "normrand", "AICzinb", "llhoodzinb", "loglikzinb"), envir = environment())
            clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b", "init.locs", "r", "var.init"), envir=environment())  
            registerDoParallel(cl)
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simNormalZINB.AIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, r, var.init)
            stopCluster(cl)

        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
            registerDoParallel(parallel::detectCores()) 
            sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simNormalZINB.AIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, r, var.init)
          
          
        } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simNormalZINB.AIC.Init(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, r, var.init)
        }
        
        cpt.loci <- sim[[1]]$loci
        
        return(list("No.BPs" = length(cpt.loci) - 2, "BP.Loc" = cpt.loci[2:(length(cpt.loci) - 1)], "AIC" = sim[[1]]$AIC, "ll" = sim[[1]]$LogLike))
      
        }
    }
}
