`GenerateSamples` <- function(n=5) {
  # generates all possible samples of size n.
Res <- NULL
for (i in 0:n) {
   AA <- i
   for (j in 0:(n-i)) {
      AB <- j
      BB <- (n-(AA+AB))
      sam <- c(AA,AB,BB)
      Res <- rbind(Res,sam)
    }
}
rownames(Res) <- 1:nrow(Res)
colnames(Res) <- c("AA","AB","BB")
return(Res)
}

