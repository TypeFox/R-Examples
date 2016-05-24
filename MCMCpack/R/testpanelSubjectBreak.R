####################################################################
## test subject-level breaks from panel residuals
## 
## written by Jong Hee Park 03/2009
## modified and integrated with other codes by JHP 07/2011
## fixed a starting.id and ending.id
######################################################################

"testpanelSubjectBreak" <-
  function(subject.id, time.id, resid, max.break=2, minimum = 10, 
           mcmc=1000, burnin=1000, thin=1, verbose=0, 
           b0, B0, c0, d0, a = NULL, b = NULL, seed = NA, 
           Time = NULL, ps.out = FALSE){   
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    N <- length(subject.id)
    
    ## groupinfo matrix
    ## col1: subj ID, col2: offset (first time C indexing), col3: #time periods
    if (min(subject.id) != 1){
      stop("subject.id should start 1!")
    }
    if (min(time.id) != 1){
      stop("time.id should start 1!")
    }
    if (is.null(Time)){
      Time <- rep(N, 1)
    }
    NC <- length(unique(subject.id))
    time.list <- as.numeric(table(subject.id))
    
    ## Make a residula list
    resid.list <- as.list(rep(NA, NC))
    start <- 1; end <- 0
    for (i in 1:NC){
      end <- start + time.list[i] - 1
      resid.list[[i]] <- ts(resid[start:end], start=Time[start])
      start <- end + 1
    }
    ## Do the break analysis
    BFout <- matrix(NA, NC, max.break + 1)
    if (ps.out ==TRUE){
      psout <- NULL
    }
    else {
      psout <- array(NA, c(max(time.list), sum(2:(max.break+1)), NC))
    }
    for (i in 1:NC){     
      residual <- resid.list[[i]]
      nk <- length(residual)
      out <- as.list(rep(NA, max.break))
      if(nk > minimum){
        for (k in 0:max.break){
          out[[k+1]] <- MCMCresidualBreakAnalysis(residual, m=k, 
                                                  b0=b0, B0=B0, c0=c0, d0=d0, a=a, b=b,
                                                  burnin=burnin, mcmc=mcmc, thin=thin, verbose=verbose, 
                                                  marginal.likelihood="Chib95")
          if (ps.out ==TRUE&k>0){
            if(k==1){
              start <- 1
            }
            else{
              start <- sum(2:k)+1
            }
            probstate <- attr(out[[k+1]], "prob.state")
            psout[1:length(probstate[,1]), start:(start+k), i] <- probstate
          }
          ## if no convergence diagnostic
          BFout[i, k+1] <- attr(out[[k+1]], "logmarglike") 
        }
      }
      if (verbose > 0){
        cat("\n ------------------------------------------------------------- ")
        cat("\n Break analysis for subject=", i, "is just finished! \n")
      }
    }
    if (ps.out ==TRUE){
      attr(BFout, "psout") <- psout
    }
    model.prob.mat <- matrix(NA, NC, max.break + 1)
    for (i in 1:NC){
      model.prob <- exp(BFout[i, ])/sum(exp(BFout[i, ]))
      winner <- which.max(model.prob)
      if (verbose > 0){
        cat("\nPr(no residual break) for subject", i, "=",
            model.prob[1])
      }
      model.prob.mat[i,] <-  model.prob
    }
    attr(BFout, "model.prob") <- model.prob.mat
    return(BFout)
  }

