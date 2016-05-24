plot.TermDocumentMatrix <- plot.DocumentTermMatrix <-
function(x,
         terms = sample(Terms(x), 20),
         corThreshold = 0.7,
         weighting = FALSE,
         attrs = list(graph = list(rankdir = "BT"),
                      node = list(shape = "rectangle", fixedsize = FALSE)),
         ...)
{
    if (system.file(package = "Rgraphviz") == "")
        stop("Plotting requires package 'Rgraphviz'.")

    m <- if (inherits(x, "TermDocumentMatrix")) t(x) else x
    m <- as.matrix(m[, terms])
    c <- cor(m)
    c[c < corThreshold] <- 0
    c[is.na(c)] <- 0
    diag(c) <- 0
    p <- Rgraphviz::plot(methods::as(c, "graphNEL"), attrs = attrs, ...)
    if (weighting) {
        i <- 1
        lw <- round(c[lower.tri(c) & c >= corThreshold] * 10)
        for (ae in Rgraphviz::AgEdge(p)) {
            Rgraphviz::lines(ae, lwd = lw[i], len = 1)
            i <- i + 1
        }
    }
    invisible(p)
}

## Plotting functions for Zipf's and Heaps'law contributed by Kurt Hornik

## See http://en.wikipedia.org/wiki/Zipf%27s_law
Zipf_plot <-
function(x, type = "l", ...)
{
    if(inherits(x, "TermDocumentMatrix"))
        x <- t(x)
    y <- log(sort(slam::col_sums(x), decreasing = TRUE))
    x <- log(seq_along(y))
    m <- lm(y ~ x)
    dots <- list(...)
    if(is.null(dots$xlab)) dots$xlab <- "log(rank)"
    if(is.null(dots$ylab)) dots$ylab <- "log(frequency)"
    do.call(plot, c(list(x, y, type = type), dots))
    abline(m)
    ## <NOTE>
    ## Perhaps this should (invisibly) return the fitted linear model
    ## instead of just the coefficients?
    coef(m)
    ## </NOTE>
}
## http://en.wikipedia.org/wiki/Heaps%27_law
## http://en.wikipedia.org/wiki/Text_corpus
## cum_vocabulary_size <-
## function(m)
## {
##     ## Should work in general, but it very slow for large simple triplet
##     ## matrices ...
##     s <- double(nrow(m))
##     v <- double(ncol(m))
##     for(i in seq_along(s)) {
##         v <- pmax(v, c(m[i, ]))
##         s[i] <- sum(v > 0)
##     }
##     s
## }
cum_vocabulary_size <-
function(m)
{
    ## Only works for simple triplet matrices.
    i <- sapply(split(m$i, m$j), min)
    tab <- table(i)
    v <- double(nrow(m))
    v[as.numeric(names(tab))] <- tab
    cumsum(v)
}

Heaps_plot <-
function(x, type = "l", ...)
{
    if(inherits(x, "TermDocumentMatrix"))
        x <- t(x)
    y <- log(cum_vocabulary_size(x))
    x <- log(cumsum(slam::row_sums(x)))
    m <- lm(y ~ x)
    dots <- list(...)
    if(is.null(dots$xlab)) dots$xlab <- "log(T)"
    if(is.null(dots$ylab)) dots$ylab <- "log(V)"
    do.call(plot, c(list(x, y, type = type), dots))
    abline(m)
    ## <NOTE>
    ## Perhaps this should (invisibly) return the fitted linear model
    ## instead of just the coefficients?
    coef(m)
    ## </NOTE>
}
