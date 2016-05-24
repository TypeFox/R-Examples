
pdols.fit <- function(data, p, icase) {

  # Data is in a list of txn matrices 
  # with the first one corresponding to y
  n <- NCOL(data[[1]])
  t <- NROW(data[[1]])
  k <- length(data) - 1
  
  if (t < 50) {
    q <- 2
  } else if (t > 50 & t <= 100) {
    q <- 3
  } else if (t > 100) {
    q <- 4
  } else {
    stop("NROW(y) must be greater than 0")
  }
  
  data1 <- lapply(data, function(z) z[(2+p):(t-p), ])
  
  ndx <- (t - 2*p - 1)
  dx <- replicate(k, matrix(0, nrow = ndx, ncol = 2*p+1), simplify = FALSE)
  
  
  for ( i in 1:n) {
    for (j in 1:(2*p + 1)) {
      for (l in 1:k) {
        dx[[l]][, j] <- data[[l+1]][(2*p+3-j):(t+1-j), i] - 
          data[[l+1]][(2*p+2-j):(t-j), i]
      }
    }
    
    for (l in 1:(k+1)) {
      data1[[l]][, i] <- data[[l]][(2+p):(t-p), i]
    }
    
    if (icase == 1) {
      xq <- do.call(cbind, dx)
    } else if (icase == 2) {
      xq <- cbind(matrix(1, nrow = ndx), do.call(cbind, dx))
    } else if (icase == 3) {
      xq <- cbind(matrix(1, nrow = ndx), c(1:ndx), do.call(cbind, dx)) 
    } else {
      stop("icase can be 1,2 or 3")
    }
    
    xq <- xq %*% solve(crossprod(xq)) %*% t(xq)
    
    for (l in 1:(k+1)) {
      data1[[l]][, i] <- data1[[l]][, i] - xq %*%  data1[[l]][, i]
    }
    
  }
  
  
  bb <- matrix(0, nrow = k, ncol = n)
  bs1 <- matrix(0, nrow = k, ncol = n)
  bs2 <- matrix(0, nrow = k, ncol = n)
  
  for (i in 1:n) {
    qqs <- do.call(cbind, lapply(data1[-1], function(z) z[, i]))
    bbs <- solve(crossprod(qqs)) %*% t(qqs) %*% data1[[1]][, i]
    rrs <- data1[[1]][, i, drop = FALSE] - qqs %*% bbs
    
    ss1 <- pdlong(rrs, q-1)
    ss2 <- qspw(rrs)
    
    
    ss1 <- ss1 * solve(crossprod(qqs)) 
    
    ss2 <- ss2 * solve(crossprod(qqs))
    ss1 <- sqrt(diag(ss1))
    ss2 <- sqrt(diag(ss2))
    
    bb[, i] <- bbs
    bs1[, i] <- ss1
    bs2[, i] <- ss2
  }
  
  v <- lapply(data1, matrix)
  vv <- do.call(cbind, v[-1])
  be <- solve(crossprod(vv)) %*% crossprod(vv, v[[1]])
  
  data2 <- lapply(data1, function(z) z - apply(t(z), 2, mean))
  
  v <- lapply(data2, matrix)
  v2 <- do.call(cbind, v[-1])
  b2 <- solve(crossprod(v2)) %*% crossprod(v2, v[[1]])
  
  
  sq1 <- sq2 <- sq3 <- sq4 <- 0
  
  for (i in 1:n) {
    www <- do.call(cbind, lapply(data1[-1], function(z) z[, i]))
    ww2 <- do.call(cbind, lapply(data2[-1], function(z) z[, i]))
    
    rr <- data1[[1]][,i] - www %*% be
    r2 <- data2[[1]][,i] - ww2 %*% b2
    
    sk1 <- pdlong(matrix(rr), q)
    sk2 <- qspw(matrix(rr))
    sk3 <- pdlong(matrix(r2), q)
    sk4 <- qspw(matrix(r2))
    
    sq1 <- sq1 + sk1 * crossprod(www)
    sq2 <- sq2 + sk2 * crossprod(www)
    sq3 <- sq3 + sk3 * crossprod(ww2)
    sq4 <- sq4 + sk4 * crossprod(ww2)
  }
  
  skk <- solve(crossprod(vv))
  skk1 <- skk %*% sq1 %*% skk
  tr01 <- sqrt(diag(skk1))
  skk2 <- skk %*% sq2 %*% skk
  tr02 <- sqrt(diag(skk2))
  
  skk <- solve(crossprod(v2))
  skk3 <- skk %*% sq3 %*% skk
  tr03 <- sqrt(diag(skk3))
  skk4 <- skk %*% sq4 %*% skk
  tr04 <- sqrt(diag(skk4))
  
  for(i in 1:k) {
    temp <- rbind(bb[i, ], bs1[i, ], bs2[i, ])
    rownames(temp) <- c(paste0("A", i), "SE.Param", "SE.Andrews")
    if (i == 1) {
      single_dols <- temp
    } else {
      single_dols <- rbind(single_dols, temp)
    }
  }
  
  pooled_ols_ordinary <- data.frame(
    PDOLS = be,
    SE.Param = tr01,
    SE.Andrews = tr02)
  
  pooled_ols_cte <- data.frame(
    PDOLS = b2,
    SE.Param = tr03,
    SE.Andrews = tr04)
  
  results <- list(`DOLS Single Equation` = single_dols,
                  `Pooled OLS CTE` = pooled_ols_cte,
                  `Pooled OLS ORD` = pooled_ols_ordinary
  )
}