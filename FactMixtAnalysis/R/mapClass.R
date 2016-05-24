mapClass <-
function (a, b) 
{
    l <- length(a)
    x <- y <- rep(NA, l)
    if (l != length(b)) {
        warning("unequal lengths")
        return(x)
    }
    aChar <- as.character(a)
    bChar <- as.character(b)
    Tab <- table(a, b)
    Ua <- dimnames(Tab)[[1]]
    Ub <- dimnames(Tab)[[2]]
    aTOb <- rep(list(Ub), length(Ua))
    names(aTOb) <- Ua
    bTOa <- rep(list(Ua), length(Ub))
    names(bTOa) <- Ub
    k <- nrow(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 1, max)
    for (i in 1:k) {
        I <- match(Max[i], Tab[i, ], nomatch = 0)
        aTOb[[i]] <- Ub[I]
    }
    if (is.numeric(b)) 
        aTOb <- lapply(aTOb, as.numeric)
    k <- ncol(Tab)
    Map <- rep(0, k)
    Max <- apply(Tab, 2, max)
    for (j in (1:k)) {
        J <- match(Max[j], Tab[, j])
        bTOa[[j]] <- Ua[J]
    }
    if (is.numeric(a)) 
        bTOa <- lapply(bTOa, as.numeric)
    list(aTOb = aTOb, bTOa = bTOa)
}
