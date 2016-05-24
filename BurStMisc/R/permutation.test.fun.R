"permutation.test.fun" <-
function (x, y = NULL, fun = function(x, y) sum(x * y), alternative = "greater", 
    trials = 1000) 
{
    if (length(y)) {
        n <- length(y)
        if (length(x) != n) 
            stop("x and y have different lengths")
        if (!is.numeric(y)) 
            stop("y must be numeric")
    }
    else {
        if (ncol(x) != 2) 
            stop("x does not have 2 columns and y is missing")
        x <- as.matrix(x)
        if (!is.numeric(x)) 
            stop("x must be numeric")
        y <- x[, 2]
        x <- x[, 1]
        n <- length(y)
    }
    if (length(alternative) != 1 || !is.character(alternative)) 
        stop("alternative must be a single character string")
    altnum <- pmatch(alternative, c("greater", "less"), nomatch = NA)
    if (is.na(altnum)) 
        stop("alternative must partially match 'greater' or 'less'")
    alternative <- c("greater", "less")[altnum]
    ranseed <- .Random.seed
    orig.score <- fun(x, y)
    if (length(orig.score) != 1) 
        stop("fun must return a single number")
    perm.scores <- numeric(trials)
    for (i in 1:trials) {
        perm.scores[i] <- fun(x, sample(y))
    }
    if (alternative == "greater") {
        extreme <- sum(perm.scores >= orig.score)
    }
    else {
        extreme <- sum(perm.scores <= orig.score)
    }
    ans <- list(original.score = orig.score, perm.scores = perm.scores, 
        stats = c(nobs = n, trials = trials, extreme = extreme,
	p.value=(extreme+1)/(trials+1)), 
        alternative = alternative, random.seed = ranseed, call = match.call())
    class(ans) <- "permtstBurSt"
    ans
}
