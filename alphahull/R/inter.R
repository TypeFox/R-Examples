inter <-
function (c11, c12, r1, c21, c22, r2) 
{
    d <- sqrt((c21 - c11)^2 + (c22 - c12)^2)
    v.x <- 0
    v.y <- 0
    theta1 <- 0
    theta2 <- 0
    if (d == 0 & r1 == r2) {
        n.cut <- Inf
    }
    else if (d > (r1 + r2)) {
        n.cut <- 0
    }
    else {
        if (d < max(r1, r2) - min(r1, r2)) {
            n.cut <- 0
        }
        else if (d == max(r1, r2) - min(r1, r2)) {
            n.cut <- 1
        }
        else if (d == r1 + r2) {
            n.cut <- 1
        }
        else {
            n.cut <- 2
            d1 <- (d^2 - r2^2 + r1^2)/(2 * d)
            d2 <- d - d1
            v.x <- (c21 - c11)/d
            v.y <- (c22 - c12)/d
            theta1 <- acos(d1/r1)
            theta2 <- acos(d2/r2)
        }
    }
    return(list(n.cut = n.cut, v1 = c(v.x, v.y), theta1 = theta1, 
        v2 = c(-v.x, -v.y), theta2 = theta2))
}
