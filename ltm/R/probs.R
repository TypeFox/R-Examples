probs <-
function (x) {
    pr <- plogis(x)
    if (any(ind <- pr == 1))
        pr[ind] <- 1 - sqrt(.Machine$double.eps)
    if (any(ind <- pr == 0))
        pr[ind] <- sqrt(.Machine$double.eps)
    pr
}
