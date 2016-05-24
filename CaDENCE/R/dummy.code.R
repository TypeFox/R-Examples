dummy.code <-
function(x) {
    if(!is.factor(x)) stop("\"x\" is not a factor")
    lvls <- levels(x)
    n.lvls <- length(lvls)
    lvls <- lvls[-n.lvls]
    categs <- matrix(0, ncol=n.lvls-1, nrow=length(x))
    for(i in seq_along(lvls)) categs[x==lvls[i],i] <- 1
    colnames(categs) <- lvls
    attr(categs, "levels") <- levels(x)
    categs
}
