lin.coef <-
function(x, y, out="both", nomiss=.8) {
  comb <- data.frame(x, y)
  comp <- subset(comb, complete.cases(comb))
  if(out=="both") {
    if(nrow(comp) >= nrow(comb)*nomiss) {
      b1 <- cor(comp[,1], comp[,2], use="pair") * (sd(comp[,2])/sd(comp[,1]))
      b0 <- mean(comp[,2]) - b1*mean(comp[,1])
    }
    else {
      b1 <- NA
      b0 <- NA
    }
  res <- cbind(b0, b1)
  }
  if(out=="slope") {
    if(nrow(comp) >= nrow(comb)*nomiss) {
      b1 <- cor(comp[,1], comp[,2], use="pair") * (sd(comp[,2])/sd(comp[,1]))
      res <- b1
    }
    else {
      b1 <- NA
    }
  res <- b1
  }
  if(out=="int") {
    if(nrow(comp) >= nrow(comb)*nomiss) {
      b1 <- cor(comp[,1], comp[,2], use="pair") * (sd(comp[,2])/sd(comp[,1]))
      b0 <- mean(comp[,2]) - b1*mean(comp[,1])
    }
    else {
      b0 <- NA
    }
  res <- b0
  }
  return(res)
}
