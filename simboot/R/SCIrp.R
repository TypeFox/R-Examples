SCIrp <-
function (x, conf.level, alternative) 
{
  DataMatrix <- x
  N <- nrow(DataMatrix)
  k <- round(conf.level * N, 0)
  RankDat <- apply(DataMatrix, 2, rank)
  switch(alternative,
         "two.sided" = {
           W1 <- apply(RankDat, 1, max)
           W2 <- N + 1 - apply(RankDat, 1, min)
           Wmat <- cbind(W1, W2)
           w <- apply(Wmat, 1, max)
           tstar <- round(sort(w)[k], 0)
           SCI <- function(x) {
             sortx <- sort(x)
             cbind(sortx[N + 1 - tstar], sortx[tstar])
           }
           SCS <- t(apply(DataMatrix, 2, SCI))
         },
         "less" = {
           W1 <- apply(RankDat, 1, max)
           tstar <- round(sort(W1)[k], 0)
           SCI <- function(x) {
             sortx <- sort(x)
             cbind(-Inf, sortx[tstar])
           }
           SCS <- t(apply(DataMatrix, 2, SCI))
         },
         "greater" = {
           W2 <- N + 1 - apply(RankDat, 1, min)
           tstar <- round(sort(W2)[k], 0)
           SCI <- function(x) {
             sortx <- sort(x)
             cbind(sortx[N + 1 - tstar], Inf)
           }
           SCS <- t(apply(DataMatrix, 2, SCI))
         })
  estimate <- apply(DataMatrix, 2, median)
  conf.int <- cbind(estimate, SCS)
  colnames(conf.int) <- cbind("estimate", "lower", "upper")
  out <- list(conf.int = conf.int, alternative = alternative)
  return(out)
}

