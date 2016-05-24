qmboots <- function(dataset, nfactors, nsteps, load="auto", rotation="varimax", indet="qindtest", fsi=TRUE, forced=T, distribution=NULL, ...) {
  startime <- Sys.time()
  nstat <- nrow(dataset)
  nqsorts <- ncol(dataset)
  # Number of iterations requested
  itnumber <- paste("Number of requested iterations:", nsteps)
  print(itnumber, quote=FALSE)
  
  #-----------------------------------------------
  # A. create objects and slots for the results
  #-----------------------------------------------
  #bootstrap results
  qmbr <- vector("list", nfactors)
  names(qmbr) <- paste0("factor_", 1:nfactors)
  n <- 1
  while (n <= nfactors) {
    qmbr[[n]] <- list()
    qmbr[[n]]$flagged <- data.frame(matrix(as.logical(NA), nrow=nqsorts, ncol=nsteps, 
                                           dimnames=list(colnames(dataset), 
                                                         paste0("step_",1:nsteps))))
    qmbr[[n]]$zsc     <- data.frame(matrix(as.numeric(NA), nrow=nstat, ncol=nsteps, 
                                           dimnames=list(row.names(dataset), 
                                                         paste0("step_",1:nsteps))))
    qmbr[[n]]$loa     <- data.frame(matrix(as.numeric(NA), nrow=nqsorts, ncol=nsteps, 
                                           dimnames=list(colnames(dataset), 
                                                         paste0("step_",1:nsteps))))
    n <- n+1
  }
  #indeterminacy tests results
  if (indet == "procrustes") {
  } 
  if (indet == "qindtest" | indet == "both") {
    #create dataframe for results and reports from INDETERMINACY issue
    qmts <- list()
    qmts[[1]] <- data.frame(matrix(NA, nrow=nfactors, ncol=nsteps, 
                                   dimnames=list(paste0("f",c(1:nfactors)), paste0("order_",1:nsteps))))
    qmts[[2]] <- data.frame(matrix(NA, nrow=nfactors, ncol=nsteps, 
                                   dimnames=list(paste0("f",c(1:nfactors)), paste0("sign_",1:nsteps))))
    names(qmts) <- c("torder", "tsign")
    qmts_log <- list()
    qmts_log[[1]] <- data.frame(matrix(NA, nrow=1, ncol=nsteps, dimnames=list("log", paste0("log_ord_",1:nsteps))))
    qmts_log[[2]] <- data.frame(matrix(NA, nrow=1, ncol=nsteps, dimnames=list("log", paste0("log_sig_",1:nsteps))))
    names(qmts_log) <- c("torder", "tsign")
  }
  
  #-----------------------------------------------
  # B. full sample Q method loadings (target matrix)
  #-----------------------------------------------
  #use manually introduced matrix of factor loadings as target, if any
  if (is.matrix(load) | is.data.frame(load)) {
    if ((nrow(load) == nqsorts) & (ncol(load) == nfactors)) {
      flagged <- qflag(loa=load, nstat=nstat)
      qm <- qzscores(dataset, nfactors, flagged=flagged, loa=load)
      target <- load
      colnames(target) <- paste0("target_f", 1:nfactors)
    } else stop("Q method input: The factor loading matrix provided in 'load' should have the correct number of Q-sort as rows and the correct number of rotated factors as columns.")
  } else if (is.character(load) & length(load) == 1) {
    if (load == "auto") {
      qm <- qmethod(dataset, nfactors, rotation=rotation, forced=forced, distribution=distribution, ...)
      flagged <- qm$flagged
      target <- as.matrix(qm$loa)
      colnames(target) <- paste0("target_f", 1:nfactors)
    } else stop("Q method input: 'load' has to be either 'auto' or a matrix.")
  }
  #create vector to build resamples
  #list of number of Q sorts, repeated as many times as bootstrap replications, to ensure that each Q-sort appears the same amount of times
  qsvector <- rep(1:nqsorts, times=nsteps)
  #generate a random order for the previous list
  qsrand <- sample(1:(nqsorts*nsteps), size=(nqsorts*nsteps), replace=FALSE)
  #order qsvector randomly
  qsvector <- qsvector[qsrand]
  
  #-----------------------------------------------
  # C. actual bootstrap
  #-----------------------------------------------
  it_count <- 1
  while (it_count <= nsteps) {
    qp <- 1 + (nqsorts*(it_count-1))
    # Create bootstrap resample
    subdata <- dataset[ , qsvector[c(qp:(qp+(nqsorts-1)))]]
    # Reshape target matrix of original factor loadings
    subtarget <-   target[qsvector[c(qp:(qp+(nqsorts-1)))],]
    # Full bootstrap step
    step_res <- qbstep(subdata=subdata, subtarget=subtarget, 
                       indet, nfactors, nqsorts, nstat, 
                       qmts=qmts, qmts_log=qmts_log, 
                       flagged=flagged, forced=forced, 
                       distribution=distribution, ...)
    # Export essential results: flagged, zsc and loa
    for (n in 1:nfactors) {
      # z-scores
      qmbr[[n]][["zsc"]][paste0("step_",it_count)] <- step_res$zsc[[n]]
      # Flagged Q-sorts and factor loadings -- important to assign to the correct columns!
      for (r in 1:nqsorts) {
        if (sum(rownames(qmbr[[n]][[3]])[r] == names(subdata)) != 0) { #if the Q sort is in the resample...
          qmbr[[n]][["flagged"]][r,it_count] <- step_res$flagged[[n]][which(rownames(qmbr[[n]][["flagged"]])[r] == names(subdata))]
          qmbr[[n]][["loa"]][r,it_count]    <- step_res$loadings[[n]][which(rownames(qmbr[[n]][["loa"]])[r] == names(subdata))]
        }
        colnames(qmbr[[n]][["flagged"]])[it_count] <- paste0("step_",it_count)
        colnames(qmbr[[n]][["loa"]])[it_count]     <- paste0("step_",it_count)
      }
    }
    # Export results of indeterminacy correction, if any
    if (indet == "both" | indet == "qindtest") {
      # Test results (logical)
      qmts[[1]][paste("order_",it_count, sep="")] <- step_res[[4]]
      qmts[[2]][paste("sign_",it_count, sep="")]  <- step_res[[5]]
      # Reports of solution implementation
      qmts_log[[1]][paste("log_ord_",it_count, sep="")] <- step_res[[6]]
      qmts_log[[2]][paste("log_sig_",it_count, sep="")] <- step_res[[7]]
    }
    # Count the iteration
    it_msg <- paste("Finished iteration number", it_count)
    print(it_msg, quote=FALSE)
    it_count <- it_count + 1
  }
  #---actual bootstrap ends here
  #-----------------------------------------------
  # D. Export indeterminacy results into one object
  #-----------------------------------------------
  qindet <- list()
  if (indet == "none") {
    qindet <- "Caution: no correction of PCA bootstrap indeterminacy issue was performed. This may introduced inflated variability in the results"
  }
  if (indet == "procrustes") {
    qindet <- "Procrustes rotation from MCMCpack was used to solve the PCA indeterminacy issue ('alignment problem')"
  }
  if (indet == "qindtest" | indet == "both") {
    qindet[[1]] <- qmts
    qindet[[2]] <- qmts_log
  }
  
  #-----------------------------------------------
  # E. Test which steps could not be swap corrected
  #-----------------------------------------------
  if (indet == "qindtest" | indet == "both") {
    errmsg <- "ERROR in ORDER swap: at least one factor in the resample is best match for two or more factors in the target"
    badsteps <- which(qindet[[2]][[1]] == errmsg)
    print(paste("Number of steps with order swap error:", 
                length(badsteps)))
    if (sum(qindet[[2]][[1]] == errmsg) == 0) badsteps <- 0
  } else badsteps <- 0
  
  #-----------------------------------------------
  # F. Summary stats for z-scores of bootstrap
  #-----------------------------------------------
  qmbs <- list() # Q method bootstrap summary
  #- - - - - - - - - - - - - - - - - - - - - - - -
  # 1. z-scores
  for (n in 1:nfactors) {
    if (sum(badsteps) > 0 & (indet == "qindtest" | indet == "both")) {
      t.zsc <- qmbr[[n]]$zsc[,-badsteps]
    } else t.zsc <- qmbr[[n]]$zsc
    qmbs[[n+1]] <- merge(describe(t(t.zsc)), 
                         t(apply(t.zsc, 1, quantile, 
                                 probs=c(0.025, 0.25, 0.75, 0.975), 
                                 na.rm=TRUE)), 
                         by='row.names', sort=FALSE)
    rownames(qmbs[[n+1]]) <- qmbs[[n+1]][,1]
    qmbs[[n+1]][,1] <- NULL
  }
  names(qmbs) <- paste0("factor", 0:nfactors)
  #- - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Factor scores (fragment adapted from qzscores.R)
  if (forced==T) qscores <- sort(dataset[,1], decreasing=FALSE)
  if (forced==F) qscores <- distribution
  
  # Build frame for factor scores
  zsc_mea <- data.frame(zsc_mea=c(1:nstat), row.names=row.names(dataset))
  zsc_bn <- data.frame(zsc_bn=c(1:nstat), row.names=row.names(dataset))
  for (f in 1:nfactors) zsc_mea[,f] <- qmbs[[f+1]]$mean
  colnames(zsc_mea) <- paste0("zsc_mea_f", c(1:nfactors))
  for (f in 1:nfactors) {
    for (s in 1:nstat) {
      # Find which statement has the current qscore rank
      statement <- order(zsc_mea[,f])[[s]]
      zsc_bn[statement,f] <- qscores[[s]]
    }
  }
  colnames(zsc_bn) <- paste0("fsc_f", c(1:nfactors))
  qmbs[[1]] <- zsc_bn
  names(qmbs)[1] <- c("Bootstraped factor scores")
  
  #-----------------------------------------------
  # G. Summary stats for Q-sort factor loadings of bootstrap
  #-----------------------------------------------
  qmbl <- list()
  for (n in 1:nfactors) {
    if (sum(badsteps) > 0 & indet == "qindtest" | indet == "both") {
      t.fla <- qmbr[[n]]$flagged[,-badsteps]
      t.loa <- qmbr[[n]]$loa[,-badsteps]
    } else {
      t.fla <- qmbr[[n]]$flagged
      t.loa <- qmbr[[n]]$loa
    }
    # Factor loadings
    qmbl[[n]] <-   merge(describe(t(t.loa)), 
                         t(apply(t.loa, 1, quantile, 
                                 probs=c(0.025, 0.25, 0.75, 0.975), 
                                 na.rm=TRUE)), 
                         by='row.names', sort=FALSE)
    rownames(qmbl[[n]]) <- qmbl[[n]][,1]
    qmbl[[n]][,1] <- NULL
    # Frequency of flagging
    qmbl[[n]]$flag_freq <- rowMeans(t.fla, na.rm=T)
    for (j in 1:nqsorts) {
      if (sum(badsteps) > 0 & indet == "qindtest" | indet == "both") {
        token <- t(qmbr[[n]][[1]][,-badsteps])[,j]
      } else token <- t(qmbr[[n]][[1]])[,j]
      qmbl[[n]][j,"flag_freq"] <- length(which(token == TRUE)) / (length(which(token == TRUE)) + length(which(token == FALSE)))
    }
    qmbl[[n]] <- qmbl[[n]][order(row.names(qmbl[[n]])), ]
    n <- n+1
  }
  
  names(qmbl) <- paste0("factor", 1:nfactors)
  #-----------------------------------------------
  # H. Export everything and report
  #-----------------------------------------------
  qmboots <- list()
  qmboots[[1]] <- qmbs
  qmboots[[2]] <- qmbr
  qmboots[[3]] <- qindet
  qmboots[[4]] <- as.data.frame(matrix(qsvector, nrow=nqsorts, 
                                       dimnames=list(NULL, paste0("bsampl_", 1:nsteps))))
  qmboots[[5]] <- qm
  qmboots[[6]] <- qscores
  qmboots[[7]] <- qmbl
  if (fsi == TRUE) {
    # Calculate FACTOR STABILITY INDEX
    fsii <- qfsi(nfactors=nfactors, nstat=nstat, qscores=qscores, zsc_bn=zsc_bn, qm=qm)
    qmboots[[8]] <- fsii
    names(qmboots) <- c("zscore-stats", "full.bts.res", "indet.tests", "resamples", "orig.res", "q.array", "loa.stats", "fsi")
  } else {
    names(qmboots) <- c("zscore-stats", "full.bts.res", "indet.tests", "resamples", "orig.res", "q.array", "loa.stats")
  }
  fintime <- Sys.time()
  duration <- paste(format(floor(difftime(fintime, startime, units="hours")[[1]]), width=2),":", format(floor(difftime(fintime, startime, units="mins")[[1]] %% 60), width=2),":", format(floor(difftime(fintime, startime, units="secs")[[1]] %% 60), width=2), sep="")
  
  cat("-----------------------------------------------\nBootstrap of ",nsteps," steps\nCall: qmboots(nfactors=",nfactors, ", nstat=", nstat, ", nqsorts=",nqsorts, ", nsteps=",nsteps, ", load=", load, ", rotation=", rotation, ", indet=",indet, ", fsi=",fsi,") \n-----------------------------------------------\nSTARTED : ",format(startime)," \nFINISHED: ", format(Sys.time()), "\nDURATION:            ", duration, " hrs:min:sec","\n-----------------------------------------------\nAlignment correction method: ",indet,"\nNumber of steps that could not be reordered: ",length(badsteps), "\n", sep="")
  
  invisible(qmboots)
}
#TO-DO: set printing method to show the bootstrapped z-scores, factor scores and stability indices