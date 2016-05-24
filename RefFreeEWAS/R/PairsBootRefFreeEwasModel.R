bootstrapPairs <- function(obs, pairID){

  pairIx <- split(obs, pairID)
  nobsTab <- table(sapply(pairIx, length))
  if(length(nobsTab)>1) {
    stop("All clusters must have the same number of observations.\n")
  }
 
  n <- length(pairIx)
  pairBoot <- pairIx[sample(1:n, n, replace=TRUE)]

  obsBoot <- unlist(lapply(pairBoot, function(u) sample(u,replace=FALSE)))
  names(obsBoot) <- NULL
  obsBoot
}

PairsBootRefFreeEwasModel <- function (mod, nboot, pairID) 
{
    BetaBoot <- array(NA, dim = c(dim(mod$Beta), 2, nboot))
    dimnames(BetaBoot)[1:2] <- dimnames(mod$Beta)
    dimnames(BetaBoot)[[3]] <- c("B", "B*")
    dimnames(BetaBoot)[[4]] <- 1:nboot
    attr(BetaBoot, "nSample") <- dim(mod$X)[1]
    for (r in 1:nboot) {
        isError <- TRUE
        while (isError) {
            catchError <- try({
                Yboot <- PairsBootOneRefFreeEwasModel(mod, pairID)
                bootFit <- RefFreeEwasModel(Yboot, mod$X, dim(mod$Lambda)[2], 
                  smallOutput = TRUE)
                BetaBoot[, , 1, r] <- bootFit$Beta
                BetaBoot[, , 2, r] <- bootFit$Bstar
            })
            isError <- inherits(catchError, "try-error")
        }
        if (r%%10 == 0) 
            cat(r, "\n")
    }
    class(BetaBoot) <- "BootRefFreeEwasModel"
    BetaBoot
}

PairsBootOneRefFreeEwasModel <- function (mod, pairID) 
{
    n2 <- dim(mod$X)[1]
    iboot <- bootstrapPairs(1:n2, pairID)
    mu <- mod$Bstar %*% t(mod$X)
    return(mu + mod$dispersion * mod$E[, iboot])
}