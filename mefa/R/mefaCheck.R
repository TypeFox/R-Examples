`mefaCheck` <-
function(x)
{
    if (!inherits(x, "mefa"))
        stop("object is not of class 'mefa'")

    mt <- rep(TRUE, 10)
    re <- list()

# when extracting or aggregating, it can be shorter, when e.g. x$samp <- NULL
# so leave it as TRUE
#    mt[1] <- length(x) == 5
    mt[2] <- is.matrix(x$xtab)
    if (!is.null(x$segm)) {
        mt[3] <- is.list(x$segm)
        if (attr(x, "nested")) {
            mt[4] <- sum(x$xtab) == sum(x$segm[[length(x$segm)]])
            } else {
            mt[4] <- sum(x$xtab) == sum(unlist(x$segm))}
        for (i in 1:dim(x)[3]) {
            mt[5] <- identical(rownames(x$xtab), rownames(x$segm[[i]]))
            mt[6] <- identical(colnames(x$xtab), colnames(x$segm[[i]]))
            }
        }
    if (!is.null(x$samp)) {
        mt[7] <- is.data.frame(x$samp)
        mt[8] <- identical(rownames(x$xtab), rownames(x$samp))
        }
    if (!is.null(x$taxa)) {
        mt[9] <- is.data.frame(x$taxa)
        mt[10] <- identical(colnames(x$xtab), rownames(x$taxa))
        }
# re list contains the final result and description of problems
    re[[1]] <- all(mt == TRUE)
    i <- 2
    if (!mt[1]) {re[[i]] <- "object length is not 5"
        i <- i + 1}
    if (!mt[2]) {re[[i]] <- "'$xtab' is not matrix"
        i <- i + 1}
    if (!mt[3]) {re[[i]] <- "'$segm' is not list"
        i <- i + 1}
    if (!mt[4]) {re[[i]] <- "sum of '$xtab' and sums in '$segm' are not equal"
        i <- i + 1}
    if (!mt[5]) {re[[i]] <- "rownames in '$xtab' and '$segm' are not identical"
        i <- i + 1}
    if (!mt[6]) {re[[i]] <- "colnames in '$xtab' and '$segm' are not identical"
        i <- i + 1}
    if (!mt[7]) {re[[i]] <- "'$samp' is not 'data.frame'"
        i <- i + 1}
    if (!mt[8]) {re[[i]] <- "rownames in '$xtab' and '$samp' are not identical"
        i <- i + 1}
    if (!mt[9]) {re[[i]] <- "'$taxa' is not 'data.frame'"
        i <- i + 1}
    if (!mt[10]) re[[i]] <- "colnames in '$xtab' and rownames in '$taxa' are not identical"

    return(re)
}

