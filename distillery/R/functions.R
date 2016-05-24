distill <- function(x, ...) {
    UseMethod("distill", x)
} # end of 'distill' function.

distill.list <- function(x, ...) {
    out <- c(unlist(x))
    return(out)
} # end of 'distill.list' function.

distill.matrix <- function(x, ...) {
    rn <- rownames(x)
    cn <- colnames(x)
    xdim <- dim(x)
    if(is.null(rn)) rn <- as.character(1:xdim[1])
    if(is.null(cn)) cn <- as.character(1:xdim[2])
    nm <- paste(rep(rn, xdim[2]), rep(cn, each=xdim[1]), sep=":")
    out <- c(x)
    names(out) <- nm
    return(out)
} # end of 'distill.matrix' function.

distill.data.frame <- function(x, ...) {
    rn <- rownames(x)
    cn <- colnames(x)
    xdim <- dim(x)
    if(is.null(rn)) rn <- as.character(1:xdim[1])
    if(is.null(cn)) cn <- as.character(1:xdim[2])
    rownames(x) <- rn
    colnames(x) <- cn
    out <- c(unlist(x))
    return(out)
} # end of 'distill.data.frame' function.

is.formula <- function(x) return(is.element("formula", class(x)))

datagrabber <- function(x, ...) {
    UseMethod("datagrabber", x)
} # end of 'datagrabber' function.

ci <- function(x, alpha=0.05, ...) {
    UseMethod("ci", x)
} # end of 'ci' function.

ci.numeric <- function(x, alpha=0.05, ...) {
    n <- length(x)
    m <- mean(x, ...)
    se.m <- sqrt(var(x, ...)/n)
    z.alpha <- qnorm(alpha/2, lower.tail=FALSE)
    out <- c(m - z.alpha * se.m, m, m + z.alpha * se.m)
    conf.level <- paste(round((1 - alpha)*100, digits=2), "%", sep="")
    names(out) <- c(paste(conf.level, " lower CI", sep=""), "mean", paste(conf.level, " upper CI", sep=""))
    attr(out, "data.name") <- deparse(substitute(x))
    attr(out, "method") <- "Normal Approximation"
    attr(out, "conf.level") <- (1 - alpha) * 100
    class(out) <- "ci"
    return(out)
} # end of 'ci.numeric' function.

ci.matrix <- function(x, alpha=0.05, ...) {
    n <- dim(x)[1]
    m <- apply(x, 2, mean, ...)
    se.m <- sqrt(apply(x, 2, var, ...)/n)
    z.alpha <- qnorm(alpha/2, lower.tail=FALSE)
    out <- cbind(m - z.alpha * se.m, m, m + z.alpha * se.m)
    conf.level <- paste(round((1 - alpha)*100, digits=2), "%", sep="")
    colnames(out) <- c(paste(conf.level, " lower CI", sep=""), "mean", paste(conf.level, " upper CI", sep=""))
    attr(out, "data.name") <- deparse(substitute(x))
    attr(out, "method") <- "Normal Approximation"
    attr(out, "conf.level") <- (1 - alpha) * 100
    class(out) <- "ci"
    return(out)
} # end of 'ci.matrix' function.

print.ci <- function(x, ...) {
    tmp <- attributes(x)
    print(tmp$data.name)
    cat("\n")
    print(tmp$method)
    if(!is.null(tmp$R)) cat(tmp$R, " iterations\n")
    cat("\n")
    if(!is.matrix(x)) {
        print(paste(names(c(x))[2], ": ", round(x[2], digits=3), sep=""))
        cat("\n")
        print(paste(tmp$conf.level, "% Confidence Interval: (", round(x[1], digits=4), ", ", round(x[3], digits=4), ")", sep=""))
    } else {
        y <- x
        attributes(y) <- NULL
        y <- matrix(y, tmp$dim[1], tmp$dim[2])
        colnames(y) <- tmp$dimnames[[2]]
        rownames(y) <- tmp$dimnames[[1]]
        print(y)
    }
    cat("\n")
    invisible()
} # end of 'print.ci' function.

