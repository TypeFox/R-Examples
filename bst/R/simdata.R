ex1data <- function(n.data, p=50)
{
    x <- runif(n.data*p)
    x <- matrix(x, ncol=p)
    p1 <- (0.97 * exp(-3*x[,1]))
    p3 <- exp(-2.5 * (x[,1]-1.2)^2)
    p2 <- (1- p1 - p3)

    # class.prob is a matrix of probabilities.
    # each column contains the probability that y is in each class for fixed x
    class.prob <- cbind(p1, p2, p3)

    u <- runif(n.data)
    y <- rep(0, n.data)

    y[u <= p1] <- 1
    y[u > p1 & u < (p1 + p2)] <- 2
    y[u >= (p1 + p2)] <- 3
    x <- as.data.frame(x)
    colnames(x) <- paste("x", 1:ncol(x), sep = "")
    return(list(x=x,y=y,p=class.prob))
}
