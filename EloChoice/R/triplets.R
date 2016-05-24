# triplets 15_06_16


triplets <- function(xdata, winner, runs=2, startvalue=0, k=100, progressbar=TRUE, mode="avg") {
  # IDs per contest
  nid <- dim(xdata)[2]
  if(max(winner) > nid) stop(paste("winner indicates", max(winner), "individuals per contest, but only", nid, "found" ))

  xdata <- apply(xdata,2,as.character)
  allids <- sort(unique(as.character(xdata)))

  if(mode=="seq") {
    # from wide version to long version with winner in first column and loser in second
    x <- do.call("rbind", sapply(1:nrow(xdata), function(X){
      cbind(xdata[X, winner[X]], sample(xdata[X, -winner[X]]))
    }, simplify=F))

    mat <- matrix(nrow = nrow(xdata), ncol = length(allids)); colnames(mat) <- allids
    ups <- dec <- numeric(nrow(mat)*(nid-1))
    curelo <- mat[1,]
    curelo[ ] <- startvalue

    i=1
    for(i in 1:nrow(x)) {
      p_win <- pnorm((curelo[x[i, 1]] - curelo[x[i, 2]])/(200 * sqrt(2)))
      (kp <- k * p_win)
      if(p_win==0.5) dec[i] <- 1
      if(p_win < 0.5) ups[i] <- 1
      # winner gets 'k-kp' points
      curelo[x[i, 1]] <- curelo[x[i, 1]] - kp + k
      curelo[x[i, 2]] <- curelo[x[i, 2]] + kp - k
      curelo; mean(curelo)
      #
    }

    res <- curelo


    if(runs>1) {
      if(progressbar) pb <- txtProgressBar(1,runs, style=3)
      for(r in 2:runs) {
        xdata <- xdata[sample(1:nrow(xdata)), ]
        x <- do.call("rbind", sapply(1:nrow(xdata), function(X){
          cbind(xdata[X, winner[X]], sample(xdata[X, -winner[X]]))
        }, simplify=F))

        ups <- dec <- numeric(nrow(mat)*(nid-1))
        curelo <- mat[1,]
        curelo[ ] <- startvalue

        for(i in 1:nrow(x)) {
          p_win <- pnorm((curelo[x[i, 1]] - curelo[x[i, 2]])/(200 * sqrt(2)))
          (kp <- k * p_win)
          if(p_win==0.5) dec[i] <- 1
          if(p_win < 0.5) ups[i] <- 1
          # winner gets 'k-kp' points
          curelo[x[i, 1]] <- curelo[x[i, 1]] - kp + k
          curelo[x[i, 2]] <- curelo[x[i, 2]] + kp - k
          curelo; mean(curelo)
          #
        }
        if(progressbar) setTxtProgressBar(pb, r)
        res <- rbind(res, curelo)
      }
      if(progressbar) close(pb)
    }

  }

  if(mode=="avg") {

    mat <- matrix(nrow = nrow(xdata), ncol = length(allids)); colnames(mat) <- allids
    ups <- dec <- numeric(nrow(mat))
    curelo <- mat[1,]
    curelo[ ] <- startvalue

    i=3
    #xdata[i, ]; xdata[i, winner[i]]; xdata[i, -winner[i]]

    for(i in 1:nrow(xdata)) {
      p_win <- pnorm((curelo[xdata[i, winner[i]]] - mean(curelo[xdata[i, -winner[i]]]))/(200 * sqrt(2)))
      (kp <- k * p_win)
      if(p_win==0.5) dec[i] <- 1
      if(p_win < 0.5) ups[i] <- 1
      # winner gets 'k-kp' points
      curelo[xdata[i, winner[i]]] <- curelo[xdata[i, winner[i]]] - kp + k
      curelo[xdata[i, -winner[i]]] <- curelo[xdata[i, -winner[i]]] + (kp - k)/ (nid-1)
      curelo; mean(curelo)
      #
    }

    res <- curelo

    if(runs>1) {
      if(progressbar) pb <- txtProgressBar(1,runs, style=3)
      for(r in 2:runs) {
        xdata <- xdata[sample(1:nrow(xdata)), ]
        ups <- dec <- numeric(nrow(mat))
        curelo <- mat[1,]
        curelo[ ] <- startvalue

        for(i in 1:nrow(xdata)) {
          p_win <- pnorm((curelo[xdata[i, winner[i]]] - mean(curelo[xdata[i, -winner[i]]]))/(200 * sqrt(2)))
          (kp <- k * p_win)
          if(p_win==0.5) dec[i] <- 1
          if(p_win < 0.5) ups[i] <- 1
          # winner gets 'k-kp' points
          curelo[xdata[i, winner[i]]] <- curelo[xdata[i, winner[i]]] - kp + k
          curelo[xdata[i, -winner[i]]] <- curelo[xdata[i, -winner[i]]] + (kp - k)/ (nid-1)
          curelo; mean(curelo)
          #
        }
        if(progressbar) setTxtProgressBar(pb, r)
        res <- rbind(res, curelo)
      }
      if(progressbar) close(pb)
    }

  }



  rownames(res) <- NULL
  return(res)
}
