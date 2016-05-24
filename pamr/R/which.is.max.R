which.is.max <- function(x)
{
        y <- seq(length(x))[x == max(x)]
        if(length(y) > 1)
                y <- sample(y, 1)
        y
}
