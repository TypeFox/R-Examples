dec2bin <-
  function(vec, npos = NA) {
    if(any(vec < 0)) 
      stop("algorithm only constructed for positive numbers.")
    binres <- NULL
    maxi <- max(vec)
    not.find.max <- TRUE
    if (all(is.na(npos))) {
      if (maxi > 0) {
        i <- 0
        while(not.find.max) {
          if((maxi / (2^i)) < 1) {
            not.find.max <- FALSE
          } else {
            i <- i + 1
          }
        } # while
        npos <- i
      } else {
        npos <- 1
      }
    }
    for(val in vec) {
      if(val > (2^npos - 1))
        stop("value cannot calculated.")
      bin <- NULL
      tmp <- val
      while(tmp != 0) {
        tmp <- val %/% 2
        bin <- c(bin, val %% 2)
        val <- tmp
      }
      if(is.null(bin)) 
        bin <- rep(0, npos)
      bin <- ifelse(bin == 0, 0, 1)
      bin2 <- rep(0, npos)
      bin2[npos:(npos - length(bin) + 1)] <- bin
      binres <- as.matrix(rbind(binres, bin2))
    } # loop for
    binres <- as.matrix(binres)
    colnames(binres) <- 2^((dim(binres)[2] - 1):0)
    rownames(binres) <- vec
    return(binres)
  }
