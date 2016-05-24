simeventhistories <- function(n, mpl, max.time, change.times = NULL, X, states.at.origin = NULL, Xstruc, 
                              partial.markov.x = NULL, partial.markov.eta = NULL){
  ## states at time origin:
  if(is.null(states.at.origin)){
    hf <- function(x, k){
      if(!is.null(x$all.to)){
        return(x[[k]])
      }
    }
    all.possible.from.states <- as.numeric(do.call(c, lapply(mpl, FUN = hf, k = 1)))
    all.first.from <- sample(all.possible.from.states, size = n, replace = TRUE)
  }else{
    all.first.from <- sample(states.at.origin, size = n, replace = TRUE)
  }
  p <- ncol(X)
  histories <- pmX <- NULL
  ## sample event histories:
  for(i in 1:n){
    history.i <- simsinglehistory(first.entry = 0,
                                  first.from = all.first.from[i], max.time, 
                                  change.times = change.times,
                                  mpl, x.i = X[i, ], 
                                  partial.markov.x, partial.markov.eta)
    pmX.i <- history.i$pmX
    history.i <- history.i$history.i
    if(!is.null(nrow(history.i))){
      history.i <- cbind(history.i, rep(i, 1))
      for(x.index in 1:p){
        history.i <- cbind(history.i, rep(X[i, x.index], 1))
      }
    }else{
      history.i <- cbind(history.i, rep(i, nrow(history.i)))
      for(x.index in 1:p){
        history.i <- cbind(history.i, rep(X[i, x.index], nrow(history.i)))
      }
    }
    histories <- rbind(histories, history.i)
    if(!is.null(pmX.i)){
      pmX <- rbind(pmX, pmX.i)
    }
    ## show progress:
    cat(".")
    if(i %% 25 == 0){
      cat(paste("  ", i, " event histories simulated.", sep=""))
      cat("\n")
    }
    rm(history.i)
  }
  ## store as msm.basics:
  trans <- paste(histories[, 3], histories[, 4], sep = "")
  msm.basics <- list("id" = histories[, 6], "entry" = histories[, 1],
                     "exit" = histories[, 2], "from" = histories[, 3],
                     "to" = histories[, 4], "trans" = trans, 
                     "delta" = histories[, 5])
  ## time-wise reduction of X:
  change.times <- c(0, change.times, Inf)
  exit <- msm.basics$exit
  X <- as.matrix(histories[, -c(1:6)])
  X.new <- NULL
  for(index in 1:ncol(Xstruc)){
    j <- unique(Xstruc[, index])
    if(length(j) > 1.5){
      x <- rep(0, nrow(X))
      for(cti in 1:(length(change.times) - 1)){
        x.sub.comp <- X[, j[cti]] * I((change.times[cti] < exit) & (exit <= change.times[cti + 1]))
        x <- x + x.sub.comp
      }
    }else{
      x <- X[, j]
    }
    X.new <- cbind(X.new, x)
  }
  ## update p:
  p <- ncol(X.new)
  ## attribute pmX in a suitable way to msm.basics and the names to Xstruc:
  if(!is.null(partial.markov.x)){
    pmX <- data.frame(pmX, row.names = 1:nrow(pmX))
    colnames(pmX) <- paste("pm.x", 1:ncol(pmX), sep = "")
    rownames(X.new) <- as.character(1:nrow(X.new))
    X.new <- cbind(X.new, pmX)
    p <- ncol(X.new)
    Xstruc <- cbind(Xstruc, matrix(nrow = nrow(Xstruc), ncol = ncol(pmX), 
                                   rep(1:ncol(pmX) + max(Xstruc), each = ncol(pmX))))
    colnames(Xstruc)[(p - ncol(pmX) + 1) : p] <- colnames(pmX)
  }
  rm(pmX)
  ## combine histories and X.new:
  for(x.index in 1:p){
    msm.basics[[7 + x.index]] <- X.new[, x.index]
    names(msm.basics)[7 + x.index] <- colnames(Xstruc)[x.index]
  }
  msm.basics <- data.frame(msm.basics)
  #print(msm.basics)
  histories <- NULL
  ## generate transition-type specific covariate versions:
  ttsce <- ttsce.names <- NULL
  for(q in sort(unique(trans))){
    for(j in 1:p){
      ttsce <- cbind(ttsce, msm.basics[, 7 + j] * I(msm.basics$trans == q))
      ttsce.names <- c(ttsce.names, paste(colnames(Xstruc)[j], q, sep = "."))
    }
  }
  ttsce <- data.frame(ttsce)
  names(ttsce) <- ttsce.names
  ntsui <- length(change.times) - 1 ## number of time sub-intervals.
  if(ntsui > 1.5){
    labels <- paste("t", 1:ntsui, sep = "")
    for(j in ttsce.names){
      for(t.index in 1:nrow(Xstruc)){
        if(t.index == 1){ ## reference category
          ttsce <- cbind(ttsce, ttsce[, j])
        }else{
          ttsce <- cbind(ttsce, ttsce[, j] * I(msm.basics$entry >= change.times[t.index]))
        }
        names(ttsce)[ncol(ttsce)] <- paste(j, rownames(Xstruc)[t.index], sep = ".")
      }
    }
  }
  ## generate transition-type indicators:
  tt.indicators <- tt.indicator.names <- NULL
  for(index in 1:length(mpl)){
    q <- c(paste(mpl[[index]]$from, mpl[[index]]$all.to, sep = ""))
    if(!is.null(mpl[[index]]$all.to)){
      for(q.index in q){
        tt.indicators <- cbind(tt.indicators, as.integer(msm.basics$trans == q.index))
        tt.indicator.names <- c(tt.indicator.names, paste("trans", q.index, sep = ""))
      }
    }
  }
  tt.indicators[which(msm.basics$delta == 0), ] <- rep(0, ncol (tt.indicators))
  tt.indicators <- data.frame(tt.indicators)
  names(tt.indicators) <- tt.indicator.names
  ## return:
  return(list(msm.basics = msm.basics, 
              ttsce = ttsce,
              tt.indicators = tt.indicators))
}