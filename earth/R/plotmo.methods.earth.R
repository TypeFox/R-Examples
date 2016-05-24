# plotmo.rpart.R: plotmo methods for earth objects

plotmo.singles.earth <- function(object, x, nresponse, trace, all1)
{
    singles <- NULL
    max.degree <- 1
    if(all1) # user wants all used predictors, not just those in degree1 terms?
        max.degree <- 99
    selected <- object$selected.terms[
                    reorder.earth(object, degree=max.degree, min.degree=1)]
    if(length(selected) > 0) {
        prednames <- object$namesx.org
        degree1.dirs <- object$dirs[selected, , drop=FALSE]
        # column numbers of dirs that have predictors in degree1 terms
        icol <- which(degree1.dirs != 0, arr.ind=TRUE)[,2]
        if(!any(sapply(x, is.factor))) # no factors in x?
            singles <- icol
        else {                         # factors in x
            colnames <- colnames(object$dirs)[icol]
            for(ipred in seq_along(prednames)) {
                if(is.factor(x[,ipred])) {
                    # This knows how to handle expanded factor names because
                    # it e.g. looks for "^pclass" in "pclass3rd"
                    # TODO this can give extra predictors if variable names alias
                    #      e.g. "x" and "x1" are both variable names
                    if(grepany(paste0("^", prednames[ipred]), colnames))
                        singles <- c(singles, ipred)
                } else if(prednames[ipred] %in% colnames)
                    singles <- c(singles, ipred)
            }
        }
        if(any(singles > length(prednames)))
            stop0("plotmo.singles.earth returned an index ",
                  "greater than the number of predictors\n",
                  "       singles=", paste(singles, collapse=","),
                  " prednames=", paste(prednames, collapse=","))
    }
    singles
}
plotmo.pairs.earth <- function(object, x, ...)
{
    pairs <- matrix(0, nrow=0, ncol=2)      # no pairs
    selected <- object$selected.terms[      # selected is all degree 2 terms
                    reorder.earth(object, degree=2, min.degree=2)]
    pairs <- vector(mode="numeric")
    for(i in selected)                      # append indices of the two preds in term i
        pairs <- c(pairs, which(object$dirs[i,] != 0))
    pairs <- unique(matrix(pairs, ncol=2, byrow=TRUE))
    if(nrow(pairs) > 0 && any(sapply(x, is.factor))) { # any columns in x are factors?
        # pairs works off expanded factor names, so replace each name
        # with index of original variable name
        # TODO this can give wrong results if variable names alias
        #      e.g. if "x" and "x1" are both variable names this takes the LAST
        #      of the matching names so correct with "x" "x1" but not "x1" "x"
        dir.colnames <- colnames(object$dirs)
        prednames <- object$namesx.org
        prednames.hat <- paste0("^", prednames)
        for(i in seq_len(nrow(pairs)))
            for(j in 1:2) {
                ipred1 <- 0
                for(ipred in seq_along(prednames.hat))
                    if(grepany(prednames.hat[ipred], dir.colnames[pairs[i, j]]))
                        ipred1 <- ipred
                if(ipred1 == 0)
                    stop0("internal error: illegal ipred1 in plotmo.pairs.earth")
                pairs[i, j] <- ipred1
            }
    }
    pairs
}
plotmo.y.earth <- function(object, trace, naked, expected.len)
{
    temp <- plotmo::plotmo.y.default(object, trace, naked, expected.len)

    # plotmo.y.default returns list(field=y, do.subset=do.subset)
    # do the same processing on y as earth does, e.g. if y is a two
    # level factor, convert it to an indicator column of 0s and 1s

    colnames <- colnames(temp$field)

    temp$field <- expand.arg(temp$field, model.env(object), trace, is.y.arg=TRUE,
                             xname=if(!is.null(colnames)) colnames else "y")

    temp
}
plotmo.pairs.bagEarth <- function(object, x, ...) # caret package
{
    pairs <- matrix(0, nrow=0, ncol=2)
    for(i in seq_along(object$fit))
        pairs <- rbind(pairs, plotmo.pairs.earth(object$fit[[i]], x))
    pairs[order(pairs[,1], pairs[,2]),]
}
plotmo.y.bagEarth <- function(object, trace, naked, expected.len)
{
    plotmo.y.earth(object, trace, naked, expected.len)
}
# back compatibility
get.plotmo.pairs.bagEarth <- function(object, env, x, trace, ...)
{
    plotmo.pairs.bagEarth(object, x, ...)
}
get.plotmo.y.bagEarth <- function(object, env, y.column, expected.len, trace)
{
    plotmo.y.bagEarth(object, trace)
}
