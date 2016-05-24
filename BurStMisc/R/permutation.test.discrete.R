"permutation.test.discrete" <-
function (x, y = NULL, scores, alternative = "greater", trials = 1000) 
{
    if (length(y)) {
        n <- length(y)
        if (length(x) != n) 
            stop("x and y have different lengths")
    }
    else {
        if (ncol(x) != 2) 
            stop("x does not have 2 columns and y is missing")
        y <- x[, 2]
        x <- x[, 1]
        n <- length(y)
    }
    x <- as.character(x)
    y <- as.character(y)
    if (length(alternative) != 1 || !is.character(alternative)) 
        stop("alternative must be a single character string")
    altnum <- pmatch(alternative, c("greater", "less"), nomatch = NA)
    if (is.na(altnum)) 
        stop("alternative must partially match 'greater' or 'less'")
    alternative <- c("greater", "less")[altnum]
    orig.tab <- table(x, y)
    otd <- dim(orig.tab)
    odnam <- dimnames(orig.tab)
    scnam <- dimnames(scores)
    if (!is.matrix(scores) || length(scnam) != 2 || !is.numeric(scores)) 
        stop("scores must be a numeric matrix with dimnames")
    scd <- dim(scores)
    if (any(scd != otd) && any(rev(scd) != otd)) {
        stop(paste("scores is not the proper size, should be", 
            otd[1], "by", otd[2]))
    }
    if (any(scd != otd)) {
        scores <- t(scores)
        scd <- dim(scores)
        scnam <- dimnames(scores)
        reverse <- TRUE
    }
    else {
        reverse <- FALSE
    }
    rownum <- match(scnam[[1]], odnam[[1]], nomatch = NA)
    if (any(is.na(rownum))) {
        if (reverse || otd[1] != otd[2]) 
            stop("bad dimnames for scores")
        scores <- t(scores)
        scd <- dim(scores)
        scnam <- dimnames(scores)
        rownum <- match(scnam[[1]], odnam[[1]], nomatch = NA)
        if (any(is.na(rownum))) 
            stop("bad dimnames for scores")
    }
    colnum <- match(scnam[[2]], odnam[[2]], nomatch = NA)
    if (any(is.na(colnum))) 
        stop("bad dimnames for scores")
    scores <- scores[rownum, colnum]
    ranseed <- .Random.seed
    orig.score <- sum(orig.tab * scores)
    perm.scores <- numeric(trials)
    for (i in 1:trials) {
        perm.scores[i] <- sum(table(x, sample(y)) * scores)
    }
    if (alternative == "greater") {
        extreme <- sum(perm.scores >= orig.score)
    }
    else {
        extreme <- sum(perm.scores <= orig.score)
    }
    ans <- list(original.score = orig.score, perm.scores = perm.scores, 
        stats = c(nobs = n, trials = trials, extreme = extreme,
	p.value=(extreme+1)/(trials+1)), alternative = alternative, 
	random.seed = ranseed, call = match.call())
    class(ans) <- "permtstBurSt"
    ans
}
