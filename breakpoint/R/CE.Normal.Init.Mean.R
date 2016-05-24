CE.Normal.Init.Mean <-
function(data, init.locs, eps=0.01, rho=0.05, M=200, h=5, a=0.8, b=0.8, distyp = 1, penalty = "mBIC", var.init = 100000, parallel=FALSE){
    
    if(is.data.frame(data) == "FALSE"| is.null(dim(data)[2])) {
      print("Error in data : dataframe only")                                                           
    } else if(dim(data)[2] != 1) {
      print("Error in data : single column dataframe only")   
    } else if(missing(init.locs)){
      print("Error: Initial locations are not provided!!!")
      
    } else { 
      
      if(distyp == 1 & penalty == "mBIC"){
        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)

        if(parallel == TRUE & .Platform$OS.type == "windows"){
          cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.sim4beta.Init.mBIC", "betarand", "fun.alpha", "fun.beta", "mBIC", "betaIntEst"), envir=environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "init.locs", "var.init"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.Init.mBIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
          stopCluster(cl)   

        } else if (parallel == TRUE & .Platform$OS.type == "unix"){ 
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.Init.mBIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
        
        } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta.Init.mBIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
        }
        
        loci.mBIC <- sim[[1]]$loci
        logLL <- llhood.MeanNormal(loci.mBIC, data, v=var(data[ ,1]), h)
        
        return(list("No.BPs" = length(loci.mBIC) - 2, "BP.Loc" = loci.mBIC[2:(length(loci.mBIC) - 1)], "mBIC value" = sim[[1]]$mBIC, "ll" = logLL))
        
      } else if(distyp == 2 & penalty == "mBIC"){
        Melite <- M*rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
          cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.simnormal.Init.mBIC", "normrand", "mBIC"), envir = environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b", "init.locs", "var.init"), envir = environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.Init.mBIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
          stopCluster(cl)   

        } else if (parallel == TRUE & .Platform$OS.type == "unix"){
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.Init.mBIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
           
        } else { 
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal.Init.mBIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
        }      

        loci.mBIC <- sim[[1]]$loci
        logLL <- llhood.MeanNormal(loci.mBIC, data, v=var(data[ ,1]), h)

        return(list("No.BPs" = length(loci.mBIC) - 2, "BP.Loc" = loci.mBIC[2:(length(loci.mBIC) - 1)], "mBIC value" = sim[[1]]$mBIC, "ll" = logLL))
        
      } else if(distyp == 1 & penalty == "BIC"){
        
        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
          cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.sim4beta.Init.Mean.BIC", "betarand", "fun.alpha", "fun.beta", "llhood.MeanNormal", "loglik.MeanNormal", "BIC.MeanNormal", "betaIntEst"), envir=environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "init.locs", "var.init"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.Init.Mean.BIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
          stopCluster(cl)   

                  } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.Init.Mean.BIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
          
          } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta.Init.Mean.BIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
        }
        
        loci.BIC <- sim[[1]]$loci
        
        return(list("No.BPs" = length(loci.BIC) - 2, "BP.Loc" = loci.BIC[2:(length(loci.BIC) - 1)], "BIC value" = sim[[1]]$BIC.Val, "ll" = sim[[1]]$LogLike))
        
      } else if(distyp == 2 & penalty == "BIC"){
        
        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
            cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.simnormal.Init.Mean.BIC", "normrand", "llhood.MeanNormal", "loglik.MeanNormal", "BIC.MeanNormal"), envir=environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b", "init.locs", "var.init"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.Init.Mean.BIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
          stopCluster(cl)   

        } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.Init.Mean.BIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
          
          } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal.Init.Mean.BIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
        }

        loci.BIC <- sim[[1]]$loci
        
        return(list("No.BPs" = length(loci.BIC) - 2, "BP.Loc" = loci.BIC[2:(length(loci.BIC) - 1)], "BIC value" = sim[[1]]$BIC.Val, "ll" = sim[[1]]$LogLike))
        
      } else if(distyp == 1 & penalty == "AIC"){

        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
          cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.sim4beta.Init.Mean.AIC", "betarand", "fun.alpha", "fun.beta", "llhood.MeanNormal", "loglik.MeanNormal", "AIC.MeanNormal", "betaIntEst"), envir=environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "init.locs", "var.init"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.Init.Mean.AIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
          stopCluster(cl)   

        } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.sim4beta.Init.Mean.AIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
          
          } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.sim4beta.Init.Mean.AIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, var.init)
        }
        
        loci.AIC <- sim[[1]]$loci
        
        return(list("No.BPs" = length(loci.AIC) - 2, "BP.Loc" = loci.AIC[2:(length(loci.AIC) - 1)], "AIC value" = sim[[1]]$AIC.Val, "ll" = sim[[1]]$LogLike))
        
      }  else if(distyp == 2 & penalty == "AIC"){

        Melite <- M * rho
        L <- length(data[, 1])
        L0 <- 1      
        k <- length(init.locs)
        
        if(parallel == TRUE & .Platform$OS.type == "windows"){
          cl <- makeCluster(parallel::detectCores(), type="SOCK") 
          clusterExport(cl, c("ce.simnormal.Init.Mean.AIC", "normrand", "llhood.MeanNormal", "loglik.MeanNormal", "BIC.MeanNormal"), envir=environment())
          clusterExport(cl, c("data", "rho", "M", "h", "eps", "Melite", "L", "L0", "a", "b", "init.locs", "var.init"), envir=environment())  
          registerDoParallel(cl)
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.Init.Mean.AIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
          stopCluster(cl)   
          
        } else if (parallel == TRUE & .Platform$OS.type == "unix"){         
          registerDoParallel(parallel::detectCores()) 
          sim <- foreach(k = k, .errorhandling = c('pass')) %dopar% ce.simnormal.Init.Mean.AIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
          
          } else {
          sim <- foreach(k = k, .errorhandling = c('pass')) %do% ce.simnormal.Init.Mean.AIC(k, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init)
        }
        
        loci.AIC <- sim[[1]]$loci
        
        return(list("No.BPs" = length(loci.AIC) - 2, "BP.Loc" = loci.AIC[2:(length(loci.AIC) - 1)], "AIC value" = sim[[1]]$AIC.Val, "ll" = sim[[1]]$LogLike))              
      }
    }
  }
