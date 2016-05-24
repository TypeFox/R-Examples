inudge.qq.plot.internal <-
function (i, obj)
{
    n <- min(which(i < apply(matrix(1:length(obj$pi)), 1, function(n) sum(obj$pi[1:n]))))
    labs <- c("unif",rep("norm",length(obj$mu)))[n];
    switch(labs, unif = runif(1, min = obj$a, max=obj$b),
      norm = rnorm(1, mean = obj$mu[n-1], sd = obj$sigma[n-1]));
}

