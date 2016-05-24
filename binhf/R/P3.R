`P3` <-
function (t) 
{
    y <- -(t^2)/(1 - t^2)
    y <- exp(y)
    y <- y * (t + 0.01)^0.25

    y <- c(y, rev(y))
	n <- length(y)
    y[seq(from=1, to=n, by=2)]

}

