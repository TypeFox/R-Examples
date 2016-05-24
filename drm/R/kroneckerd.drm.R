"kroneckerd.drm" <-
  function (y, nclass = length(levels(factor(y))), dep = "I") 
{
  nrep <- length(y)
  mis <- which(is.na(y) == TRUE)
  if (length(mis) > 0) {
    if (any(!is.na(y[min(mis):nrep]))) 
      stop("Only monotone missing patterns allowed when imposing a selection model")
    y <- y[seq(min(mis))]
    mis <- min(mis)
    ymis <- as.matrix(expand.grid(rep(list(1:nclass), length(mis))))
    y <- rep(list(y), nrow(ymis))
    for (i in 1:length(y)) {
      y[[i]][mis] <- ymis[i, ]
    }
    nre <- length(y[[1]])
    if (length(grep("M", dep)) > 0) {
      if (length(grep("M2", dep)) > 0) {
        ma <- matrix(0, ncol = nrep - max(nre, 3), nrow = nclass^3)
        ma[1, ] <- 1
        if (nre == 2) 
          wmiss <- sapply(y, function(k, nclass) kronecker.drm(c(k, 
                                                                NA), nclass), nclass = nclass)
        else wmiss <- sapply(y, function(k, nclass, nre, 
                                         nrep) {
          cbind(sapply(1:(nre - 2), function(j, k, nclass) kroneckerd.drm(k[j:(j + 
                                                                              2)], nclass = nclass), nclass = nclass, k = k), 
                ma, if (nrep > 3) {
                  if (nre < 4) 
                    ma
                  else cbind(sapply(2:(nre - 2), function(j, 
                                                          k, nclass) c(kronecker.drm(k[j:(j + 1)], 
                                                                                    nclass = nclass), rep(0, nclass^3 - nclass^2)), 
                                    nclass = nclass, k = k), ma)
                })
        }, nclass = nclass, nre = nre, nrep = nrep)
        wmiss <- array(wmiss, dim = c(nclass^3, if (nrep > 
                                2) 2 * nrep - 5 else 1, nclass^(length(mis))))
        wmiss
      }
      else {
        ma <- matrix(0, ncol = nrep - max(nre, 2), nrow = nclass^2)
        ma[1, ] <- 1
        wmiss <- sapply(y, function(k, nclass, nre, nrep) {
          cbind(sapply(1:(nre - 1), function(j, k, nclass) kroneckerd.drm(k[j:(j + 
                                                                              1)], nclass = nclass), nclass = nclass, k = k), 
                ma, if (nrep > 2) {
                  if (nre < 3) 
                    ma
                  else cbind(sapply(2:(nre - 1), function(j, 
                                                          k, nclass) c(kronecker.drm(k[j], nclass = nclass), 
                                                                       rep(0, nclass^2 - nclass)), nclass = nclass, 
                                    k = k), ma)
                })
        }, nclass = nclass, nre = nre, nrep = nrep)
        wmiss <- array(wmiss, dim = c(nclass^2, if (nrep > 
                                2) 2 * nrep - 3 else 1, nclass^(length(mis))))
        wmiss
      }
    }
    else {
      wmiss <- sapply(y, function(k, nclass, nrep, nre) c(kroneckerd.drm(k, 
                                                                        nclass = nclass), rep(0, (nclass^nrep) - (nclass^nre))), 
                      nclass = nclass, nrep = nrep, nre = nre)
      wmiss
    }
  }
  else {
    ymat <- matrix(0, ncol = nrep, nrow = nclass)
    ymat[, y == 1] <- c(1, rep(-1, (nclass - 1)))
    for (i in 2:nclass) {
      ymat[, y == i] <- c(rep(0, (i - 1)), 1, rep(0, (nclass - 
               i)))
    }
    if (length(grep("M", dep)) > 0) {
      wmiss <- kronecker.drm(y, nclass = nclass, dep = dep)
      array(wmiss, dim = c(dim(wmiss), 1))
    }
    else {
      w <- eval(parse(text = paste(paste(rep("c(", nrep - 
                        1), collapse = ""), paste("ymat[,1]", paste("%*% t(ymat[,", 
                                                                    2:nrep, "]))", collapse = "")))))
      matrix(w)
    }
  }
}

