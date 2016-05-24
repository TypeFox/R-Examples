# Conversion from BaBooN's BBPMM to coda's mcmc or mcmc.list object, or to mice's mids object (almost)
# Version:            0.1
# Date:        2013-12-17
# Author: T.S., ctb: F.M.
# Note:   Needs Hmisc's asNumericMatrix and matrix2dataframe
#         and coda's mcmc and mcmc.list function.
# Further infos, references and credits:
# See for mice: Van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate 
#               Imputation by Chained Equations in R. Journal of Statistical Software,
#               Vol. 45, No. 3, pp. 1--67. URL http://www.jstatsoft.org/v45/i03/.
#       and for the mids-class the help page to "mids" / "mids-class" / "Multiply imputed data set" within mice. 
#
# See for coda: Plummer, M. and Best,N. and Cowles, K. and Vines, K. (2006) CODA: Convergence 
#               Diagnosis and Output Analysis for MCMC, R News, Vol. 6, pp. 7--11
#       and for the mcmc-class the help page to "mcmc", for mcmc.list the help page of
#       "mcmc.list" within coda.
#
# See for Hmisc: Harrell, F.E., with contributions from Charles Dupont and many others. (2013) Hmisc: Harrell
#                Miscellaneous. R package version 3.13-0. http://CRAN.R-project.org/package=Hmisc
# License: GPL-2 | GPL-3 | (GPL >= 2)

# impdiagnosticconversion ----------------------------------------------------
# Version:            0.1
# Date:        2013-12-17
# Author: T.S., ctb: F.M.
# Note:   Needs Hmisc's asNumericMatrix and matrix2dataframe
#         and coda's mcmc and mcmc.list function.
# License: GPL-2 | GPL-3 | (GPL >= 2)


impdiagnosticconversion <- function(imputed.data,
                            type=c("mcmc.list","mcmc","mids")){

  if(class(imputed.data) != "imp" || imputed.data$call[[1]] != "BBPMM"){
    stop("Only BaBooN's BBPMM-function atm.")
  } 
  
  type <- match.arg(type)
  
  ## Conversion of chains
  m                    <- imputed.data$M
  var_names            <- names(imputed.data$Chains)
  num_imps             <- length(imputed.data$Chains)
  chainconversion_out  <- impChainConversion(imputed.data$Chains,
                                 m,
                                 var_names,
                                 num_imps)
    
  if(type != "mids"){
    
    clen1 <- dim(chainconversion_out[[1]])[3]
    clen2 <- dim(chainconversion_out[[2]])[3]
    clen3 <- dim(chainconversion_out[[3]])[3]
    clen4 <- dim(chainconversion_out[[4]])[3]
    
    coda_out       <- list()
    coda_out$means <- vector("list",clen1)
    coda_out$vars  <- vector("list",clen2)
    coda_out$medians <- vector("list",clen3)
    coda_out$sds     <- vector("list",clen4)
    
    for(i in 1:clen1){
      coda_out$means[[i]] <- mcmc(t(chainconversion_out$mean_arr[,,i]))
    }
    
    for(i in 1:clen2){
      coda_out$vars[[i]]  <- mcmc(t(chainconversion_out$var_arr[,,i]))
    }
    
    for(i in 1:clen3){
      coda_out$medians[[i]] <- mcmc(t(chainconversion_out$median_arr[,,i]))
    }
    
    for(i in 1:clen4){
      coda_out$sds[[i]]  <- mcmc(t(chainconversion_out$sd_arr[,,i]))
    }
    
    
    if(type == "mcmc.list"){
      
      coda_out$means   <- mcmc.list(coda_out$means)
      coda_out$vars    <- mcmc.list(coda_out$vars)
      coda_out$medians <- mcmc.list(coda_out$medians)
      coda_out$sds     <- mcmc.list(coda_out$sds)

    }  

    return(coda_out)
    
  } else {
  
  
  mids_out        <- vector("list",length=17)
  
  names(mids_out) <- c("call", "data",
                       "m", "nmis", "imp",
                       "method", "predictorMatrix",
                       "visitSequence", "form",
                       "post", "seed", "iteration", "lastSeedValue",
                       "chainMean", "chainVar",
                       "loggedEvents", "pad")

  ### Conversion of (meta) information
  mids_out$call          <- imputed.data$call
  mids_out$nmis          <- imputed.data$mis.num
  mids_out$m             <- imputed.data$M
  mids_out$seed          <- ifelse(is.null(imputed.data$seed),NA, imputed.data$seed)
  mids_out$iteration     <- imputed.data$nIter
  mids_out$lastSeedValue <- imputed.data$LastSeed
  mids_out$method        <- rep("pmm",length(imputed.data$Chains))

  ### BBPMM provides less information atm, hence empty list elements
  mids_out$pad             <- list()
  mids_out$predictorMatrix <- NULL
  mids_out$post            <- character(0)
  mids_out$visitSequence   <- integer(0)
  mids_out$loggedEvents    <- NULL
  
  ### Retrieving chains
  mids_out$chainMean <- chainconversion_out$mean_arr
  mids_out$chainVar  <- chainconversion_out$var_arr
  
  ### START data and imp conversion

  ## data
  mids_out$data <- asNumericMatrix(imputed.data$impdata[[1]])
  is.na(mids_out$data[,var_names][which(imputed.data$indMatrix[,var_names])]) <- TRUE
  
  mids_out$data <- matrix2dataFrame(mids_out$data)

  ## imp
  
  all_var_names <- names(imputed.data$mis.num)
  help_list <- vector("list",length(all_var_names))
  
  for(i in 1:length(all_var_names)){
    if(imputed.data$mis.num[i] == 0){
      help_list[[i]] <- NULL
      } else {
        help_list[[i]]            <- sapply(imputed.data$impdata,function(x) {x[[i]]})
        help_list[[i]]            <- as.matrix(help_list[[i]][imputed.data$indMatrix[,i],])
        row.names(help_list[[i]]) <- which(imputed.data$indMatrix[,i])
        colnames(help_list[[i]])  <- 1:m
      }
  }
  
  names(help_list) <- all_var_names
  mids_out$imp     <- help_list
  
  ### END data and imp conversion
  
  ### visitSequence... 
  
  mids_out$visitSequence <- vector("integer",num_imps)
  
  for(j in 1:num_imps){
    if(any(var_names[j]==all_var_names)){
      mids_out$visitSequence[j] <- which(var_names[j]==all_var_names)
    }
  }
  
  names(mids_out$visitSequence) <- var_names


  ### Change class to "mids"

  class(mids_out) <- "mids"
  return(mids_out)
  }

}

# ####################################################################### #
# impChainConversion -------------------------------------------------
# Auxiliary function for impdiagnosticconversion
# Necessary, because BBPMM stores all imputed values. The means, variances
# medians and standard deviations get calculated here.
# Hence, it is expandable if other statistics are needed apart
# from these four.
# Prepares chains
# Version:           0.1
# Date:       2013-12-17
# Author:           T.S.
# License: GPL-2 | GPL-3 | (GPL >= 2)

impChainConversion <- function(chai, m, var_names,num_imps){
  
    mv_list              <- vector("list", 4)
    names(mv_list)       <- c("means", "variance","median", "sd")
    mv_list              <- lapply(mv_list, function(x) x <- vector("list", 
                                   length = num_imps))
    names(mv_list$means) <-
      names(mv_list$variance) <-
      names(mv_list$median) <- 
      names(mv_list$sd) <- var_names
    
    
    for (i in 1:num_imps) {
        mv_list$means[[i]]    <- lapply(chai[[i]], colMeans)
        mv_list$variance[[i]] <- lapply(chai[[i]], function(x) apply(x, 
                                        2, var))
        mv_list$median[[i]] <- lapply(chai[[i]], function(x) apply(x, 
                                        2, median))
        
        mv_list$sd[[i]] <- lapply(chai[[i]], function(x) apply(x, 
                                        2, sd))
        
        #expand here for different statistics calculations
    }

    iterLens <- sapply(mv_list$means[[1]], length)
    if (sum(abs(diff(iterLens))) != 0) {
        stop("Conversion stopped. Chains have different lengths.")
    }
    
    help_list          <- vector("list", 4)
    help_list$mean     <- help_list$var <- help_list$median <- help_list$sd <- list()
    mean_arr           <- array(0, dim = c(num_imps, iterLens[1], m))
    var_arr            <- array(0, dim = c(num_imps, iterLens[1], m))
    median_arr         <- array(0, dim = c(num_imps, iterLens[1], m))
    sd_arr             <- array(0, dim = c(num_imps, iterLens[1], m))
    
    dimnames(mean_arr) <-
      dimnames(var_arr) <-
      dimnames(median_arr) <-
      dimnames(sd_arr) <-
      list(var_names, 
           1:iterLens[1], paste("Chain", 1:m))
    
    #..and expand the loops for different statistics
    for (j in 1:m) {
        help_list$mean[[j]] <- 
           help_list$var[[j]] <-
           help_list$median[[j]] <- 
           help_list$sd[[j]] <- matrix(0,
                                       num_imps,
                                       iterLens[1])
        
        for (k in 1:num_imps) {
            help_list$mean[[j]][k, ]   <- mv_list$means[[k]][[j]]
            help_list$var[[j]][k, ]    <- mv_list$variance[[k]][[j]]
            help_list$median[[j]][k, ] <- mv_list$median[[k]][[j]]
            help_list$sd[[j]][k, ]     <- mv_list$sd[[k]][[j]]
        }
        mean_arr[, , j] <- help_list$mean[[j]]
        var_arr[, , j]  <- help_list$var[[j]]
        median_arr[, , j]  <- help_list$median[[j]]
        sd_arr[, , j]   <- help_list$sd[[j]]
    }
    
    return(list(mean_arr = mean_arr,
                var_arr  = var_arr,
                median_arr = median_arr,
                sd_arr     = sd_arr))

}
