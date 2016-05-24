multExpFDR <-function (theScores, thePermScores, type = "pos") 
{
    numberOfPermutations <- length(thePermScores)
    if (!type %in% c("two.sided", "pos", "neg")) 
        stop("Wrong type!")
    ff <- function(x) -x
    biggerEq <- function(x, y) {
        y <- sort(y, decreasing = TRUE)
        a <- match(x, x)
        b <- x %in% y
        d <- match(x, sort(c(x, y), decreasing = TRUE))
        return(d - a + b)
    }
    theFDR <- matrix(NA, nrow = nrow(theScores), ncol = ncol(theScores))
    rownames(theFDR) <- rownames(theScores)
    if (type == "two.sided") {
        theScores <- abs(theScores)
        thePermScores <- lapply(thePermScores, abs)
    }
    if (type == "neg") {
        theScores <- -theScores
        thePermScores <- lapply(thePermScores, ff)
    }
    for (i in 1:ncol(theScores)) {
        ord <- order(theScores[, i], decreasing = TRUE)
        z <- theScores[ord, i]
        randomZ <- as.vector(sapply(thePermScores, function(x) x[, 
            i]))
        randomZ <- sort(randomZ, decreasing = TRUE)
        numberisBigger <- biggerEq(z, randomZ)
        theFDR[ord, i] <- numberisBigger/((1:length(z)) * numberOfPermutations)
    }
    return(theFDR)
}