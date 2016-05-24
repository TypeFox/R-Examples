rDist <- function(x, ...) {
    d <- density(na.omit(x), ...)
    
    structure(function(n)
        rnorm(n = n,
              mean = sample(x = x, size = n, replace = TRUE),
              sd = d$bw),
              class = "rDist")
}
plot.rDist <- function(x, ...) {
    eval(expression(plot(d, ...)), envir = environment(x))
}
