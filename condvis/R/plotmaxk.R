plotmaxk <- 
function (maxk)
{
    q <- quantile(maxk, probs = seq(0, 1, 0.1))
    plot(seq(0, 1, 0.1), q, type = "l", xlab = "proportion of data", ylab = 
        "max k attained")
    points(seq(0, 1, 0.1), q, pch = 16)
}

update.maxk <- 
function (object, ...)
{

}