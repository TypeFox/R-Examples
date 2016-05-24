
 library(glmbb)
 data(crabs)

 gout <- glmbb(satell ~ color * spine * width * weight,
     data = crabs, cutoff = Inf)

 fits <- ls(envir = gout$envir, pattern = "^sha1")

 form1 <- Map(function (x) get(x, envir = gout$envir)$formula, fits)

 form2 <- sapply(form1, tidy.formula.hierarchical)
 names(form2) <- NULL

 foo <- Vectorize(function(x) isHierarchical(as.formula(x)))
 all(foo(form2))

 library(digest)

 digest.interaction.graph <- function(formula) {
    stopifnot(inherits(formula, "formula"))
    mt <- terms(formula)
    mf <- attr(mt, "factors")
    if (! is.matrix(mf)) {
        g <- "empty graph"
    } else {
        mr <- attr(mt, "response")
        if (mr != 0)
            mf <- mf[- mr, , drop = FALSE]
        g <- matrix(0, nrow(mf), nrow(mf))
        for (i in 1:nrow(mf))
            for (j in 1:nrow(mf))
                if (i != j)
                    g[i, j] <- any(mf[i, ] * mf[j, ] >= 1)
        dimnames(g) <- list(rownames(mf), rownames(mf))
        idx <- match(sort(rownames(mf)), rownames(mf))
        g <- g[idx, , drop = FALSE]
        g <- g[ , idx, drop = FALSE]
    }
    digest(g, algo = "sha1")
 }

 nterms <- function(formula) {
    stopifnot(inherits(formula, "formula"))
    mt <- terms(formula)
    ml <- attr(mt, "term.labels")
    length(ml)
 }

 foo <- Vectorize(function(x) isGraphical(as.formula(x)))
 is.graphical.1 <- foo(form2)

 foo <- Vectorize(function(x) digest.interaction.graph(as.formula(x)))
 graphs <- foo(form2)
 fred <- split(form2, graphs)
 foo <- Vectorize(function(x) nterms(as.formula(x)))
 fred.max <- sapply(fred, function(x) x[foo(x) == max(foo(x))])
 is.graphical.2 <- form2 %in% fred.max
 all(is.graphical.1 == is.graphical.2)

