`stcs` <-
function(dframe, expand = FALSE, drop.zero = FALSE, zero.pseudo="zero.pseudo")
{
    x <- dframe
    if (any(is.na(x)))
        stop("'xtab' contains 'NA'")
    if (min(dim(as.matrix(x))) == 1)
        stop("2-4 columns required")
    nrows <- nrow(x)
    ncols <- ncol(x)
    if (ncols > 4) stop("2-4 columns required")
    if (length(zero.pseudo) > 2)
        stop("length of 'zero.pseudo' should be <= 2")
    if (length(zero.pseudo) == 1)
        zero.pseudo <- rep(zero.pseudo, 2)
# evaluate according to input column numbers
    x <- data.frame(x)
    if (ncols == 2)
        x <- data.frame(as.factor(x[, 1]), as.factor(x[, 2]),
            rep(1, nrows), rep("undefined", nrows))
    if (ncols == 3) {
        if (is.factor(x[,3]) || is.character(x[,3])) {
            x <- data.frame(x[,1:2], rep(1, nrows), x[,3])
            } else {
            x <- data.frame(x, rep("undefined", nrows))}}
    if (!is.numeric(x[, 3]))
        stop("count column must be numeric")
# treating zero cases
    if (any(x[, 3] == 0)) {
        zpart <- x[x[, 3] == 0,]
        x <- x[x[, 3] != 0,]
        } else {
        zpart <- NULL}
    if (drop.zero && !is.null(zpart))
        zpart <- NULL
# expand argument
    if (expand) {
#        tmp <- data.frame(inflate(x[, c(1:2,4)], x[, 3])) ## inflate is obsolete
        tmp <- rep.data.frame(x[, c(1:2,4)], x[, 3])
        x <- data.frame(tmp[, 1:2], rep(1, sum(x[, 3])), tmp[, 3])}
# check match with predefined characters
    if (!is.null(zpart)){
        if (any(zero.pseudo %in% unique(as.character(x[, 2]))))
            stop("'zero.pseudo' found in taxa names: specify other value")
        if (any(zero.pseudo %in% unique(as.character(x[, 4]))))
            stop("'zero.pseudo' found in segment names: specify other value")
        if ("not.defined" %in% unique(as.character(x[, 2])))
            stop("'not.defined' found in taxa names: change the name")
        if ("not.defined" %in% unique(as.character(x[, 4])))
            stop("'not.defined' found in segment names: change the name")
        if ("not.defined" %in% zero.pseudo)
            stop("'not.defined' found in 'zero.pseudo': change the argument")
        zpart[,2] <- rep(zero.pseudo[1], nrow(zpart))
        zpart[,4] <- rep(zero.pseudo[2], nrow(zpart))
        if (identical(zero.pseudo[1], zero.pseudo[2]))
            zero.pseudo <- zero.pseudo[1]
# gives warning for mismatch of zero.pseudo and count=0 cases
        joint <- which(zpart[,1] %in% x[,1])
        if (length(joint) != 0) {
            zpart <- zpart[-joint,]
            if (nrow(zpart) == 0)
                zpart <- NULL
            warning("zero count for non zero sample found")}
        colnames(x) <- 1:4
        colnames(zpart) <- 1:4
# put the 2 parts together
        x <- merge(x, zpart, all = TRUE)
        } else zero.pseudo <- "not.defined"
# refresh dimnames
    rownames(x) <- NULL
    colnames(x) <- c("samp", "taxa", "count", "segm")
# set type
    x$samp <- as.factor(x$samp)
    x$taxa <- as.factor(x$taxa)
    x$segm <- as.factor(x$segm)
# drop unused levels
    x[] <- lapply(x, function(x) x[drop = TRUE])
    if (all(x$count %in% c(0, 1)))
        expand <- TRUE
    class(x) <- c("stcs", "data.frame")
    attr(x, "call") <- match.call()
    attr(x, "expand") <- expand
    attr(x, "zero.count") <- any(x$count == 0)
    attr(x, "zero.pseudo") <- zero.pseudo
    return(x)
}
