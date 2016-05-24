twoSampleLinearRankTest <-
function (x, y, location.shift.null = 0, scale.shift.null = 1, 
    alternative = "two.sided", test = "wilcoxon", shift.type = "location") 
{
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    test <- match.arg(test, c("wilcoxon", "normal.scores", "moods.median", 
        "savage.scores"))
    shift.type <- match.arg(shift.type, c("location", "scale"))
    data.name <- c(deparse(substitute(x)), deparse(substitute(y)))
    names(data.name) <- c("x", "y")
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if ((bad.obs.x <- sum(!(ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[ok]
        warning(paste(bad.obs.x, "observations with NA/NaN/Inf in 'x' removed."))
    }
    if (length(unique(x)) < 2) 
        stop("'x' must contain at least two distinct non-missing observations.")
    if (!is.vector(y, mode = "numeric") || is.factor(y)) 
        stop("'y' must be a numeric vector")
    if ((bad.obs.y <- sum(!(ok <- is.finite(y)))) > 0) {
        is.not.finite.warning(y)
        y <- y[ok]
        warning(paste(bad.obs.y, "observations with NA/NaN/Inf in 'y' removed."))
    }
    if (length(unique(y)) < 2) 
        stop("'y' must contain at least two distinct non-missing observations.")
    if (shift.type == "location" && !missing(location.shift.null)) {
        if ((length(location.shift.null) != 1) || !is.finite(location.shift.null)) 
            stop("'location.shift.null' must be a single finite numeric value")
        x <- x - location.shift.null
    }
    if (shift.type == "scale" && !missing(scale.shift.null)) {
        if ((length(scale.shift.null) != 1) || !is.finite(scale.shift.null) || 
            scale.shift.null <= 0) 
            stop("'scale.shift.null' must be a single finite positive numeric value")
        x <- x/scale.shift.null
    }
    m <- length(x)
    n <- length(y)
    N <- m + n
    R <- 1:N
    scores <- switch(test, wilcoxon = (2/(N + 1)) * R - 1, normal.scores = qnorm(R/(N + 
        1)), moods.median = sign(R - (N + 1)/2), savage.scores = cumsum(1/(N - 
        R + 1)))
    z <- c(x, y)
    index <- c(rep(1, m), rep(2, n))
    index <- index[order(z)]
    z <- sort(z)
    R <- rank(z)
    if (any(table(R) > 1)) {
        vec <- sapply(split(scores, R), mean)
        scores <- vec[match(R, as.numeric(names(vec)))]
    }
    V <- sum(scores[index == 1])
    E.V <- m * mean(scores)
    Var.V <- ((m * n)/N) * var(scores)
    z <- (V - E.V)/sqrt(Var.V)
    p.value <- switch(alternative, two.sided = 2 * (1 - pnorm(abs(z))), 
        less = pnorm(z), greater = 1 - pnorm(z))
    stat <- z
    names(stat) <- "z"
    parameters <- NULL
    string <- switch(alternative, two.sided = "!=", less = "<", 
        greater = ">")
    null.value <- "Fx(t)"
    if (shift.type == "location") 
        names(null.value) <- ifelse(location.shift.null == 0, 
            "Fy(t)", paste("Fy(t - ", location.shift.null, ")", 
                sep = ""))
    else {
        names(null.value) <- ifelse(scale.shift.null == 1, "Fy(t)", 
            paste("Fy(t / ", scale.shift.null, ")", sep = ""))
    }
    alternative <- paste(names(null.value), string, "Fx(t) for at least one t")
    string1 <- switch(test, wilcoxon = "Wilcoxon Rank Sum Test", 
        normal.scores = "Normal Scores Test", savage.scores = "Savage Scores Test", 
        moods.median = "Mood's Median Test")
    method <- paste("Two-Sample Linear Rank Test:", string1, 
        "Based on Normal Approximation", sep = paste("\n", space(33), 
            sep = ""))
    ret.list <- list(statistic = stat, parameters = parameters, 
        p.value = p.value, estimate = NULL, null.value = null.value, 
        alternative = alternative, method = method, estimation.method = NULL, 
        sample.size = c(nx = m, ny = n), data.name = data.name, 
        bad.obs = c(x = bad.obs.x, y = bad.obs.y))
    oldClass(ret.list) <- "htest"
    ret.list
}
