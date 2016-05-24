gng.qq.plot.internal <-
function (i, obj)
{
    n <- min(which(i < apply(matrix(1:length(obj$pi)), 1, function(n) sum(obj$pi[1:n]))))
    labs <- c("negexp", rep("norm", length(obj$mu)), "posexp")[n]
    switch(labs, negexp = -(rgamma(1, shape = 1, scale = obj$b[1]) +
        obj$th1), posexp = rgamma(1, shape = 1, scale = obj$b[2]) +
        obj$th2, norm = rnorm(1, obj$mu[n - 1], obj$sigma[n -
        1]))
}

