points.w.dups <-
function (x, ..., method = c("jitter", "number", "standard"), 
    factor.x = 1, factor.y = 1, dup.cex = 0.85 * par("cex")) 
{
    method <- match.arg(method)
    if (method == "standard") {
        points(x, ...)
    }
    else {
        ldots <- list(...)
        nam.ldots <- names(ldots)
        len.ldots <- length(ldots)
        if (!is.vector(x) || len.ldots == 0) 
            stop(paste("When 'method' is not 'standard'", "you must supply 'x' and 'y' vectors"))
        if (len.ldots == 1 || is.null(nam.ldots) || nam.ldots[1] == 
            "y") 
            index <- 1
        else {
            index <- match("y", nam.ldots)
            if ((is.na(index))) 
                index <- 1
        }
        y <- ldots[[index]]
        if (!is.vector(y) || length(y) != length(x)) 
            stop(paste("'y' must be the second argument after 'x', or else", 
                "you must name it.  Also, 'y' must be a vector that", 
                "is the same length as 'x'"))
        ldots <- ldots[-index]
        tab <- table(x, y)
        if (!any(tab > 1)) {
            arg.list <- c(list(x = x, y = y), ldots)
            do.call("points", arg.list)
        }
        else {
            d <- dimnames(tab)
            x.vals <- as.numeric(d[[1]])
            y.vals <- as.numeric(d[[2]])
            x.mat <- matrix(x.vals, nrow = nrow(tab), ncol = ncol(tab))
            y.mat <- matrix(y.vals, nrow = nrow(tab), ncol = ncol(tab), 
                byrow = TRUE)
            if (any(index <- tab == 1)) {
                arg.list <- c(list(x = x.mat[index], y = y.mat[index]), 
                  ldots)
                do.call("points", arg.list)
            }
            if (method == "jitter") {
                amount.x <- factor.x * min(diff(sort(unique(x))), 
                  na.rm = TRUE)/5
                amount.y <- factor.y * min(diff(sort(unique(y))), 
                  na.rm = TRUE)/5
                index <- (as.character(x) %in% x.mat[tab > 1]) & 
                  (as.character(y) %in% y.mat[tab > 1])
                dum.x <- x[index]
                dum.y <- y[index]
                arg.list <- c(list(x = jitter(dum.x, amount = amount.x), 
                  y = jitter(dum.y, amount = amount.y)), ldots)
                do.call("points", arg.list)
            }
            else {
                for (i in unique(tab[tab > 1])) {
                  arg.list <- c(list(x = x.mat[tab == i], y = y.mat[tab == 
                    i], labels = as.character(i), cex = dup.cex), 
                    ldots)
                  do.call("text", arg.list)
                }
            }
        }
    }
}
