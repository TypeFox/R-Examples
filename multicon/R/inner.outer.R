inner.outer <-
function(L) {
  DF <- data.frame(L[1])
  for(i in 2:length(L)) {DF <- cbind(DF, data.frame(L[i]))}
  MAT <- cor(DF, use="pair")
  cols <- unlist(lapply(L, ncol))
  startindx <- c(0,cumsum(cols)) + 1
  endindx <- startindx - 1
  INNER <- rep(NA, length(L))
  OUTER <- rep(NA, length(L))
  for(i in 1:length(L)) {
    INNER[i] <- mean(MAT[startindx[i]:endindx[i+1],startindx[i]:endindx[i+1]][upper.tri(MAT[startindx[i]:endindx[i+1],startindx[i]:endindx[i+1]])])
    OUTER[i] <- mean(MAT[-startindx[i]:-endindx[i+1],startindx[i]:endindx[i+1]])
  }
  res <- rbind(INNER, OUTER)
  res.names <- rep(NA, length(L))
  for(i in 1:length(L)) {res.names[i] <- paste("Set", i, sep="")}
  colnames(res) <- res.names
  rownames(res) <- c("Inner r", "Outer r")
  return(res)
}
