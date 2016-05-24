
isa.in.silico <- function(num.rows=300, num.cols=50, num.fact=3,
                          mod.row.size=round(.5*num.rows/num.fact),
                          mod.col.size=round(.5*num.cols/num.fact),
                          noise=0.1,
                          mod.signal=rep(1, num.fact),
                          mod.noise=rep(0, num.fact),
                          overlap.row=0, overlap.col=overlap.row) {

  isa.status("Creating in-silico data set", "in")
  
  if (max(mod.row.size) > num.rows || max(mod.col.size) > num.cols) {
    stop("Inconsistent data configuration")
  }

  if (length(mod.noise) != num.fact) {
    stop("Invalid `mod.noise' length")
  }
  
  mod.row.size <- rep(mod.row.size, length=num.fact)
  mod.col.size <- rep(mod.col.size, length=num.fact)

  ol.row <- overlap.row * seq(0,num.fact-1)
  ol.col <- overlap.col * seq(0,num.fact-1)
  crow.from <- cumsum(c(1,mod.row.size))[-(num.fact+1)] - ol.row
  crow.to   <- cumsum(mod.row.size) - ol.row
  ccol.from <- cumsum(c(1,mod.col.size))[-(num.fact+1)] - ol.col
  ccol.to   <- cumsum(mod.col.size) - ol.col

  mod.signal <- rep(mod.signal, length=num.fact)
  
  rowMod <- matrix(0, nrow=num.rows, ncol=num.fact)
  colMod <- matrix(0, nrow=num.cols, ncol=num.fact)

  data <- matrix(0, nrow=num.cols, ncol=num.rows)

  for (i in seq_len(num.fact)) {
    rowMod[crow.from[i]:crow.to[i],i] <- sqrt(mod.signal[i])
    colMod[ccol.from[i]:ccol.to[i],i] <- sqrt(mod.signal[i])
  }
  
  data[] <- colMod %*% t(rowMod)
  data[] <- data + rnorm(length(data), mean=0, sd=noise)
  rowMod[] <- ifelse(rowMod != 0, 1, 0)
  colMod[] <- ifelse(colMod != 0, 1, 0)

  ## Add module specific noise
  for (i in which(mod.noise != 0)) {
    idx2 <- crow.from[i]:crow.to[i]
    idx1 <- ccol.from[i]:ccol.to[i]
    rnd <- rnorm(length(idx1)*length(idx2), sd=mod.noise[i])
    data[idx1,idx2] <- data[idx1,idx2] + rnd
  }
  
  res <- list(data=t(data), rowMod=rowMod, colMod=colMod)

  isa.status("DONE", "out")

  res
}

ppa.in.silico <- function(num.rows1=300, num.rows2=200, num.cols=50,
                          num.fact=3,
                          mod.row1.size=round(.5*num.rows1/num.fact),
                          mod.row2.size=round(.5*num.rows2/num.fact),
                          mod.col.size=round(.5*num.cols/num.fact),
                          noise=0.1,
                          mod.signal=rep(1, num.fact),
                          mod.noise=rep(0, num.fact),
                          overlap.row1=0, overlap.row2=overlap.row1,
                          overlap.col=overlap.row1) {

  isa.status("Creating in-silico data set for PPA", "in")

  if (max(mod.row1.size) > num.rows1 ||
      max(mod.row2.size) > num.rows2 ||
      max(mod.col.size) > num.cols) {
    stop("Inconsistent data configuration")
  }

  if (length(mod.signal) != num.fact) {
    stop("Invalid `mod.signal' length")
  }
  if (length(mod.noise) != num.fact) {
    stop("Invalid `mod.noise' length")
  }

  mod.row1.size <- rep(mod.row1.size, length=num.fact)
  mod.row2.size <- rep(mod.row2.size, length=num.fact)
  mod.col.size  <- rep(mod.col.size,  length=num.fact)

  ol.row1 <- overlap.row1 * seq(0, num.fact-1)
  ol.row2 <- overlap.row2 * seq(0, num.fact-1)
  ol.col  <- overlap.col  * seq(0, num.fact-1)
  crow1.from <- cumsum(c(1,mod.row1.size))[-(num.fact+1)] - ol.row1
  crow1.to   <- cumsum(mod.row1.size)-ol.row1
  crow2.from <- cumsum(c(1,mod.row2.size))[-(num.fact+1)] - ol.row2
  crow2.to   <- cumsum(mod.row2.size)-ol.row2
  ccol.from <- cumsum(c(1,mod.col.size))[-(num.fact+1)] - ol.col
  ccol.to   <- cumsum(mod.col.size)-ol.col

  mod.signal <- rep(mod.signal, length=num.fact)

  row1Mod <- matrix(0, nrow=num.rows1, ncol=num.fact)
  row2Mod <- matrix(0, nrow=num.rows2, ncol=num.fact)
  colMod <- matrix(0, nrow=num.cols, ncol=num.fact)

  data1 <- matrix(0, nrow=num.cols, ncol=num.rows1)
  data2 <- matrix(0, nrow=num.cols, ncol=num.rows2)

  for (i in seq_len(num.fact)) {
    row1Mod[crow1.from[i]:crow1.to[i],i] <- sqrt(mod.signal[i])
    row2Mod[crow2.from[i]:crow2.to[i],i] <- sqrt(mod.signal[i])
    colMod[ccol.from[i]:ccol.to[i],i] <- sqrt(mod.signal[i])
  }

  data1[] <- colMod %*% t(row1Mod)
  data1[] <- data1[] + rnorm(length(data1), mean=0, sd=noise)  
  data2[] <- colMod %*% t(row2Mod)
  data2[] <- data2[] + rnorm(length(data2), mean=0, sd=noise)
  row1Mod[] <- ifelse(row1Mod != 0, 1, 0)
  row2Mod[] <- ifelse(row2Mod != 0, 1, 0)
  colMod[] <- ifelse(colMod != 0, 1, 0)

  ## Add module specific noise
  for (i in which(mod.noise != 0)) {
    idx2 <- crow1.from[i]:crow1.to[i]
    idx3 <- crow2.from[i]:crow2.to[i]
    idx1 <- ccol.from[i]:ccol.to[i]
    rnd <- rnorm(length(idx1)*length(idx2), sd=mod.noise[i])
    data1[idx1,idx2] <- data1[idx1,idx2] + rnd
    rnd <- rnorm(length(idx1)*length(idx3), sd=mod.noise[i])
    data2[idx1,idx3] <- data2[idx1,idx2] + rnd
  }

  res <- list(data1=t(data1), data2=t(data2),
              row1Mod=row1Mod, row2Mod=row2Mod, colMod=colMod)

  isa.status("DONE", "out")

  res
}
