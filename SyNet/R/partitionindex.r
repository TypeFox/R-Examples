partitionindex <- function (iptsymp, replica = 1)
{
    if (!is.matrix(iptsymp)) {
        cat("Error: the input argument is not a matrix \n")
        return(invisible())
    }
    if (!isSymmetric(iptsymp)) {
        cat("Error: input matrix is not symmetric \n")
        return(invisible())
    }
    replica <- as.integer(replica)
    if (replica < 1 | replica > 5000) {
        cat("Error: the argument 'replica' is outside the interval [1,5000] \n")
        return(invisible())
    }
    iptsymp <- ifelse(iptsymp > 0, 1, 0)
    diag(iptsymp) <- 0
    n <- nrow(iptsymp)
    np <- choose(n, 2)
    probtie <- mean(iptsymp[row(iptsymp) > col(iptsymp)])
    Partition_Index <- function (msym) {
      coef <- array(0, n)
      cperf <- array(Inf, n)
      for (i in 1:n) {
          neigh <- which(msym[i, ] == 1)
          a <- length(neigh)
          if (a > 1) coef[i] <- sum(msym[neigh, neigh])/(a * (a - 1))
          cperf[neigh] <- pmin(cperf[neigh], coef[i]) #In this loop I update
                                                      #the clustering performance
                                                      #at the neighborhood
      }
      Idx <- sum(pmax(coef, cperf))/n
      return(Idx)
    }
    realval <- Partition_Index(iptsymp)
    if (replica == 1) {
        cat("Partition Index of Sympatry Matrix: ", realval, "\n")
        return(invisible(realval))
    }
    rndval <- c()
    i <- 1
    a <- 0
    rnd <- matrix(0, n, n)
    while (i <= replica) {
        aux <- ifelse(runif(np) <= probtie, 1, 0)
        rnd[] <- 0
        rnd[row(rnd) > col(rnd)] <- aux
        rnd <- rnd + t(rnd)
        rndval <- c(rndval, Partition_Index(rnd))
        if (rndval[i] >= realval) a <- a + 1
        i <- i + 1
    }
    hist(rndval, main = paste("RANDOM TEST \n", "N: ", replica,
        "; PI (p >= PI): ", round(realval, 3), "(", round(a/replica,
            3), ")"), xlab = "Partition Index (PI)", col = 1)
    abline(v = realval, col = 2)
    diag(iptsymp) <- 1
    return(list(ProbTie = probtie, PIobserved = realval,
           PIrandomized = fivenum(rndval), TestValue = round(a/replica, 5)))
}
