arrange.jm <- function (x, par=FALSE) 
{
  marks <- names(x)
  #marks <- gsub("S1_", "*Scaff", marks)
  dis <- dim(x)
  rows <- dis[1]
  cols <- dis[2]
  newm <- matrix(0, nrow = cols, ncol = rows + 2)
  newm[1:cols, 3:(rows + 2)] <- t(x)
  newm[which(is.na(newm))] <- "-"
  for (i in 1:cols) {
    if (length(which(newm[i, c(1:30)] %in% c("nn", "np") == TRUE)) > 0) {
      newm[i, 2] = "<nnxnp>;"
    }
    if (length(which(newm[i, c(1:30)] %in% c("lm", "ll") == TRUE)) > 0) {
      newm[i, 2] = "<lmxll>;"
    }
    if (length(which(newm[i, c(1:30)] %in% c("hh", "hk", "kk") == TRUE)) > 0) {
      newm[i, 2] = "<hkxhk>;"
    }
    if (length(which(newm[i, c(1:30)] %in% c("ab", "cd", "ac", "ad", "bc", "bd") == TRUE)) > 0) {
      newm[i, 2] = "<abxcd>;"
    }
    #if (newm[i, c(1:3)] == "-") {
    #  newm[i, 2] = "-"
    #}
  }
  newm[, 1] <- marks
  #newm[1:10,1:10]
  #v <- which(newm[, 2] == "-")
  #z <- which(newm[, 2] == 0)
  #z <- unique(c(v, w))
  if(par==F){
    newm <- newm[, -c(3, 4)]
  } else{newm <- newm}
  newm[, 3:(dim(newm)[2] - 1)] <- paste(newm[, 3:(dim(newm)[2] - 
                                                    1)], sep = "")
  rownames(newm) <- newm[, 1]
  newm3 <- as.data.frame(newm[, -1])
  namesjoin <- matrix(paste(rownames(newm3[, ]), newm3[, 1]), 
                      nrow = length(rownames(newm3)), ncol = 1)
  calls <- newm3[, 2:(dim(newm3)[2])]
  dims <- dim(calls)
  joinmap <- matrix(NA, nrow = dims[1] * 2, ncol = dims[2])
  v <- seq(1, dims[1])
  w <- seq(2, 2 * dims[1], by = 2)
  q <- seq(1, (2 * dims[1]), by = 2)
  for (i in 1:dims[1]) {
    v2 <- v[i]
    w2 <- w[i]
    q2 <- q[i]
    joinmap[w2, ] <- as.matrix(calls[v2, ])
    joinmap[q2, 1] <- as.matrix(namesjoin[v2, 1])
  }
  return(joinmap)
}