perks.d <-
function(x, a, b, p, xsc){
    x <- x - xsc
    f1 <- 1 + exp(b - p*x)
    mod <- a/f1
    db <- (1/f1-1)*a/f1
    dp <- - db*x
    .aArg <- as.list(match.call()['a'])
    if (unlist(lapply(.aArg, is.name))){
        .grad <- array(0, c(length(x), 3L), list(NULL, c('a', 'b', 'p')))
        da <- 1/f1
        .grad[,'a'] <- da
    } else .grad <- array(0, c(length(x), 2L), list(NULL, c('b', 'p')))
    .grad[,'b'] <- db
    .grad[,'p'] <- dp
    attr(mod, "gradient") <- .grad
    attr(mod, 'xsc') <- xsc
    mod
}
