`upd.prob` <-
function(samples, v, weight, comp.list)
{
    p <- matrix(0, nrow=nrow(v), ncol=ncol(v))
    s <- nrow(samples)
    p <- apply(samples,2,function(x) 
        table(x)[match(as.character(comp.list), dimnames(table(x))[[1]])]/s)
    p[is.na(p)] <- 0
    (1-weight)*v + weight*p
}

