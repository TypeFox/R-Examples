stratify <-
function (dataset){
    xx <- apply(dataset, 1, function(x) paste(x, collapse = "\r"))
    tab <- table(xx)
    st <- names(tab)
    strata <- match(xx,st)
    return(strata)
}
