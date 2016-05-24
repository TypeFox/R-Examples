table.qtukey <- function(alpha, n)
{
    u <- matrix(nrow = n, ncol = 4)
    colnames(u) <- c("alpha/2", "alpha", "1-alpha", "1-alpha/2")
    for(i in 2:n)
    {
        a <- function(i) {
            qtukey(alpha / 2, i, Inf)
        }
        b <- function(i) {
            qtukey(alpha, i, Inf)
        }
        g <- function(i) {
            qtukey(1 - alpha, i, Inf)
        }
        d <- function(i) {
            qtukey(1 - (alpha / 2), i, Inf)
        }
        u[i, ] <- c(a(i), b(i), g(i), d(i))
    }
    y <- u[2:n, ]
    rownames(y) <- c(2:n)
    print(y)
}