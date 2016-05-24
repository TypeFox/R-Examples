pseudo.v <- function (data, statistic, order = 1, matrix = FALSE) 
{
    if (matrix == TRUE) 
        n <- nrow(data)
    if (matrix == FALSE) 
        n <- nrow(as.matrix(data))
    theta.star <- statistic(data)
    theta.star.i <- matrix(ncol = 1, nrow = n)
    Pseudo.v <- matrix(ncol = 1, nrow = n)
    for (i in 1:n) {
        if (matrix == TRUE) {
		if(order == 1) theta.star.i[i] <- statistic(data[-i, ])
		if(order > 1) theta.star.i[i] <- statistic(data[-c(i:(order + i - 1)),])
        }
        if (matrix == FALSE) {
            data <- as.matrix(data)
            if(order == 1) theta.star.i[i] <- statistic(data[-i])
		if(order > 1) theta.star.i[i] <- statistic(data[-c(i:(order + i - 1))])
        }
        Pseudo.v <- n * theta.star - ((n - 1) * theta.star.i)
    }
    data.frame(theta.hat.star = theta.star.i, Pseudo.val = Pseudo.v)
}