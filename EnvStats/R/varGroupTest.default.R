varGroupTest.default <-
function (object, group, test = "Levene", correct = TRUE, data.name = NULL, 
    group.name = NULL, parent.of.data = NULL, subset.expression = NULL, 
    ...) 
{
    test <- match.arg(test, c("Levene", "Bartlett"))
    if (is.null(data.name)) 
        data.name <- deparse(substitute(object))
    if (is.null(group.name)) 
        group.name <- deparse(substitute(group))
    x <- as.vector(unlist(object))
    if (!is.numeric(x)) 
        stop("All elements of 'object' must be numeric")
    if (!is.factor(group)) 
        group <- factor(unlist(group))
    if (any(!is.finite(group))) 
        stop("NA's/Inf's not allowed in 'group'.")
    if (length(x) != length(group)) 
        stop("'object' and 'group' must have the same number of elements.")
    x.list <- split(x, group)
    n.grps <- length(x.list)
    if (n.grps == 1) {
        warning("Only one group supplied, so the function 'varTest' was called.")
        ret.list <- varTest(x = unlist(x.list), data.name = data.name)
    }
    else {
        c.names <- names(x.list)
        if (identical(c.names, as.character(1:n.grps))) {
            c.names <- paste(group.name, c.names, sep = ".")
        }
        bad.obs.vec <- numeric(n.grps)
        if (na.inf.flag <- any(!is.finite(x))) {
            for (i in 1:n.grps) {
                x.i <- x.list[[i]]
                if ((bad.obs <- sum(!(x.ok <- is.finite(x.i)))) > 
                  0) {
                  is.not.finite.warning(x.i)
                  x.i <- x.i[x.ok]
                  warning(paste(bad.obs, "observations with NA/NaN/Inf in the group", 
                    c.names[i], "removed."))
                  if (length(x.i) < 2 || all(x.i == x.i[1])) 
                    stop(paste("There must be at least 2", "unique non-missing values in each group. ", 
                      "This is not true for the following group: ", 
                      c.names[i]))
                  bad.obs.vec[i] <- bad.obs
                  x.list[[i]] <- x.i
                }
            }
        }
        names(bad.obs.vec) <- c.names
        sample.size.vec <- sapply(x.list, length)
        names(sample.size.vec) <- c.names
        estimate <- sapply(x.list, var)
        names(estimate) <- c.names
        sep.string <- paste("\n", space(33), sep = "")
        switch(test, Levene = {
            if (na.inf.flag) {
                x <- unlist(x.list)
                group <- factor(rep(c.names, sample.size.vec))
            }
            z <- abs(lm(x ~ factor(group))$residuals)
            aov.list <- anova(lm(z ~ group))
            statistic <- aov.list[1, "F value"]
            names(statistic) <- "F"
            parameters <- aov.list[, "Df"]
            names(parameters) <- c("num df", "denom df")
            p.value <- aov.list[1, "Pr(>F)"]
            method <- paste("Levene's Test for", "Homogenity of Variance", 
                sep = sep.string)
        }, Bartlett = {
            df <- n.grps - 1
            nu.vec <- sample.size.vec - 1
            nu <- sum(nu.vec)
            sigma.sq <- sum((nu.vec/nu) * estimate)
            statistic <- nu * log(sigma.sq) - sum(nu.vec * log(estimate))
            if (correct) {
                C <- 1 + (1/(3 * df)) * sum(1/nu.vec - 1/nu)
                statistic <- statistic/C
            }
            names(statistic) <- "Chisq"
            parameters <- c(df = df)
            p.value <- 1 - pchisq(statistic, df = df)
            string <- paste(ifelse(correct, "(With", "(Without"), 
                "Correction Factor)")
            method <- paste("Bartlett's Test for", "Homogenity of Variance", 
                string, sep = sep.string)
        })
        null.value <- 1
        names(null.value) <- "Ratio of each pair of variances"
        ret.list <- list(statistic = statistic, parameters = parameters, 
            p.value = p.value, estimate = estimate, null.value = null.value, 
            alternative = paste("At least one variance differs"), 
            method = method, data.name = data.name, grouping.variable = group.name)
        if (!is.null(parent.of.data)) 
            ret.list$parent.of.data <- parent.of.data
        if (!is.null(subset.expression)) 
            ret.list$subset.expression <- subset.expression
        ret.list <- c(ret.list, list(sample.size = sample.size.vec, 
            bad.obs = bad.obs.vec))
        oldClass(ret.list) <- "htest"
    }
    ret.list
}
