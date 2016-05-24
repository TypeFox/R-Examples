backShift <- function(X, ExpInd, covariance=TRUE, ev=0, threshold =0.75, 
                      nsim=100, sampleSettings=1/sqrt(2), 
                      sampleObservations=1/sqrt(2), nodewise=TRUE, 
                      tolerance = 10^(-4), baseSettingEnv = 1, verbose = FALSE){

    # checking the validity of the input arguments
    if(!is.matrix(X) & !is.data.frame(X)) 
      stop("'X' must be a matrix of data frame")

    X <- as.matrix(X)
    p <- ncol(X)
    
    if(!is.vector(ExpInd)) 
      stop("ExpInd must be a vector")
    if( length(ExpInd) != nrow(X)) 
      stop("'ExpInd' must have as many entries as X has rows")
    if(ev<0) 
      stop("ev must be non-negative")
    if( threshold <=0.5 | threshold >1) 
      stop("threshold must be between 0.5 and 1")
    if(nsim < 2) 
      stop("'nsim' must be at least 2 (but usually at least 20-100")
    if(length((settings <- unique(ExpInd)))<2) 
      stop("need at least three different settings")
     if(sampleSettings<=0) 
       stop("sampleSettings needs to be positive")
    if(sampleObservations<=0) 
      stop("sampleObservations needs to be positive")
    if(sampleSettings>1) 
      stop("sampleSettings needs to be at most 1")
    if(sampleObservations>1) 
      stop("sampleObservations needs to be at most 1")
    
    # stability selection parameters
    q <- if(!nodewise) sqrt(ev*(2*threshold-1)*(p^2-p)) else sqrt(ev*(2*threshold-1))
    subs <- sampleSettings* length(settings)
    drawE <- function(x){
        z <- floor(x)
        if( runif(1) <  x-z) z <- z+1
        z <- max(3,z)
        return(z)
    }        
    
    ## point estimator
    
    # run FFDIAG
    Deltalist <- computeDelta(X,  ExpInd, covariance=covariance)$Delta
    Delta <- array(unlist(Deltalist), 
                   dim = c(nrow(Deltalist[[1]]), 
                           ncol(Deltalist[[1]]), length(Deltalist)))
    tryCatch({
      estimatedB <- ffdiag(Delta,eps=tolerance, itermax = 500)$B  
    }, 
    warning=function(w){
      warning("WARNING: backShift -- point estimate -- diagonalization ", 
              "did not succeed -- result not trustworthy", immediate. = TRUE)
    },
    error=function(e){
      warning("ERROR: backShift -- point estimate -- diagonalization did not ", 
              "succeed. Possible model mispecification, returning the empty graph.", 
              immediate. = TRUE)
    },
    finally = {
      estimatedB <- try(suppressWarnings(
        ffdiag(Delta,eps=tolerance, itermax = 500)$B), silent = TRUE)
      if(inherits(estimatedB, "try-error"))
        return(list(Ahat=0*diag(p), 
                    AhatAdjacency = 0*diag(p), 
                    varianceEnv = matrix(0, nrow = length(unique(ExpInd)), ncol = p)))
    }
    )
    
    # permute and scale
    res.point.est <- permuteAndScale(estimatedB, verbose)
    Ahat <- res.point.est$Ahat
    
    # estimate intervention variances
    varEnv <- computeVarEnv(res.point.est$rescaledDhat, 
                            Deltalist, baseSettingEnv, verbose)

    # if ev > 0: run stability selection
    if(ev>0){
        if(verbose){
          cat("Starting stability selection... \n")
        }
      
        AhatList <- vector("list", nsim)
        numberOfRunsNotConverged <- 0
        
        for(i in 1:nsim){
          if(verbose){
            cat("Stability selection: Iteration", i, "... \n")
          }
          
          useSettings <- sample(settings, drawE(subs))
          useSamples <- NULL
          for(s in 1:length(useSettings)){
            ind <- which(ExpInd %in% useSettings[s])
            useSamples <- 
              c(useSamples, 
                sort(sample(ind, round(length(ind)*sampleObservations))))
          }
      
          Xcurrent <- X[useSamples,]
          ExpIndcurrent <- ExpInd[useSamples]
          
          ## difference between the covariance matrices of X
          ## under the different interventions
          Deltalist <- computeDelta(Xcurrent, 
                                    ExpIndcurrent, 
                                    covariance=covariance)$Delta
          Delta <- array(unlist(Deltalist), 
                         dim = c(nrow(Deltalist[[1]]), 
                                 ncol(Deltalist[[1]]), length(Deltalist)))
                    
          rm(estimatedB)
          tryCatch({
            estimatedB <- ffdiag(Delta,eps=tolerance, itermax = 500)$B
          }, 
            warning=function(w){
              numberOfRunsNotConverged <<- numberOfRunsNotConverged + 1
            },
            error=function(e){
              warning("ERROR: backShift -- stability selection -- ", 
                      "diagonalization did not succeed. Possible model ", 
                      "mispecification or sample size too small. ",
                      "Returning the point estimates only. The adjacency matrix is empty.",
                      immediate. = TRUE)
            },
            finally = {
              estimatedB <- 
                try(suppressWarnings(ffdiag(Delta,eps=tolerance, itermax = 500)$B), 
                    silent = TRUE)
              if(inherits(estimatedB, "try-error"))
                return(list(Ahat=Ahat, 
                            AhatAdjacency = 0*diag(p), 
                            varianceEnv = varEnv))
            }
          )     
                  
          Ahat.stab.sel <- permuteAndScale(estimatedB, verbose)$Ahat            
          AhatList[[i]] <- edgeSelection(Ahat.stab.sel, q,  nodewise=nodewise)
        }

        if(verbose){
          cat("Stability selection...done! \n")
        }
        
        numberOfRunsConverged <- 100*(nsim-numberOfRunsNotConverged)/nsim
        
        message("backShift: Percentage of runs in stability selection that converged: ", 
                numberOfRunsConverged, "%")
        if(numberOfRunsConverged < 75)  
          warning("WARNING: backShift -- stability selection -- ",
                  "only ", numberOfRunsConverged,"% of the runs converged", 
                  immediate. = TRUE)
        
        AhatAdjacency <- edgeRetention(AhatList, threshold, p)
    }else{
        AhatAdjacency <- NULL
    }
    list(Ahat=Ahat, AhatAdjacency = AhatAdjacency, varianceEnv = varEnv)
}
