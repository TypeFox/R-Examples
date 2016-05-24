pamr.knnimpute.old <- function(data, k = 10) {
  x <- data$x
  N <- dim(x)
  p <- N[2]
  
  N <- N[1]
        col.nas  <- apply(x, 2, is.na)
  if ((sum(col.nas) == N) > 0) {
    stop("Error: A column has all missing values!")
  }
  nas <- is.na(drop(x %*% rep(1, p)))
  xcomplete <- x[!nas,  ]
  xbad <- x[nas,,drop=FALSE ]
  xnas <- is.na(xbad)
  xbadhat <- xbad
  cat(nrow(xbad), fill = TRUE)
  for(i in seq(nrow(xbad))) {
    cat(i, fill = TRUE)
    xinas <- xnas[i,  ]
    xbadhat[i,  ] <- nnmiss(xcomplete, xbad[i,  ], xinas, K = k)
  }
  x[nas,  ] <- xbadhat
  data2 <-data
  data2$x <-x
  return(data2)
}

nnmiss <- function(x, xmiss, ismiss, K = 1) {
  xd <- scale(x, xmiss, FALSE)[, !ismiss]
  dd <- drop(xd^2 %*% rep(1, ncol(xd)))
  od <- order(dd)[seq(K)]
  xmiss[ismiss] <- drop(rep(1/K, K) %*% x[od, ismiss, drop = FALSE])
  xmiss
}

