`balanced.specaccum` <-
function (comm, permutations=100, strata=strata, grouped=TRUE, reps=0, scale=NULL) {
    accumulator <- function(x, ind) {
        rowSums(apply(x[ind, ], 2, cumsum) > 0)
    }
    stratified.sample <- function(factor,grouped=TRUE,reps=0) {
        n <- length(factor)
        levs <- levels(droplevels(factor))
        minimum <- min(summary(factor))
        if (reps > 0) {
            alllevs <- summary(factor)
            goodlevs <- alllevs > (reps-1)
            levs <- names(alllevs[goodlevs])
            minimum <- reps
        }
        nl <- length(levs)
        seq2 <- array(nl*minimum)
        seq1 <- sample(n)
        strat <- sample(nl)
        count <- 0
        for (i in 1:nl) {
            for (j in 1:n) {
                if (factor[seq1[j]]==levs[strat[i]]) {
                        count <- count+1
                        if (count > i*minimum) {count <- count-1} 
                        seq2[count] <- seq1[j]
                }
            }
        }
        if (grouped==FALSE) {
            seq3 <- sample(seq2)
            seq2 <- seq3
        }
        return(seq2)
    }
    x <- comm
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (p == 1) {
        x <- t(x)
        n <- nrow(x)
        p <- ncol(x)
    }
    specaccum <- sdaccum <- sites <- perm <- NULL
    if (n == 1) 
        stop(paste("only 1 site provided"))
    if (is.factor(strata) != TRUE) 
        stop(paste("strata should be a categorical variable"))
    n1 <- length(stratified.sample(strata,grouped,reps))
    perm <- array(dim = c(n1, permutations))
    for (i in 1:permutations) {        
         perm[, i] <- accumulator(x, stratified.sample(strata,grouped,reps))
    }
    sites <- 1:n1
    specaccum <- apply(perm, 1, mean)
    sdaccum <- apply(perm, 1, sd)
    out <- list(call = match.call(), method = "balanced species accumulation", sites = sites, 
        richness = specaccum, sd = sdaccum, perm = perm)
    class(out) <- "specaccum"
    if (is.null(scale)!=TRUE) {
        n <- length(strata)
        levs <- levels(droplevels(strata))
        if (reps > 0) {
            alllevs <- summary(strata)
            goodlevs <- alllevs > (reps-1)
            levs <- names(alllevs[goodlevs])
        }
        nlevs <- length(levs)
        tot <- 0
        for (i in 1:nlevs) {
            ind <- strata==levs[i]
            tot <- tot + mean(scale[ind])
        }
        tot <- tot/nlevs
        out$sites <- out$sites * tot
    }
    out
}

