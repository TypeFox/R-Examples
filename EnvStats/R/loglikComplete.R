loglikComplete <-
function (theta, x, distribution = "norm") 
{
    sum(log(do.call(paste("d", distribution, sep = ""), args = c(list(x = x), 
        as.list(theta)))))
}
