`vark` <-
function(x, y)
    {
        ties.x <- rle(sort(x))$lengths
        ties.y <- rle(sort(y))$lengths
        n <- length(x)
        t1 <- n * (n - 1) * (2 * n + 5)
        t2 <- sum(ties.x * (ties.x - 1) * (2 * ties.x + 5))
        t3 <- sum(ties.y * (ties.y - 1) * (2 * ties.y + 5))
        v1 <- (t1 - t2 - t3)/18
        t1 <- sum(ties.x * (ties.x - 1) * (ties.x - 2))
        t2 <- sum(ties.y * (ties.y - 1) * (ties.y - 2))
        v2 <- (t1 * t2)/(9 * n * (n - 1) * (n - 2))
        t1 <- sum(ties.x * (ties.x - 1)) * sum(ties.y * (ties.y - 1))
        v3 <- t1/(2 * n * (n - 1))
        v1 + v2 + v3
    }

