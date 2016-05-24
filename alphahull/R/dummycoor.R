dummycoor <-
function (tri.obj, l1, l2, m, away) 
{
    v <- l2 - l1
    v <- c(v[2], -v[1])
    norm <- sum(v^2)
    if (norm > 0) {
        v <- v/norm
    }
    mp <- (l1 + l2)/2
    eps <- 1e-05
    test <- mp + eps * v
    inconv <- in.convex.hull(tri.obj, test[1], test[2])
    if (inconv) {
        dum <- m - away * v
    }
    else {
        dum <- m + away * v
    }
    return(dum)
}
