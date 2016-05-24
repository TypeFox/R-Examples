require(fracdiff)

## confirm that we guessed right:
## fracdiff.sim(....., d = 0, backComp = FALSE)  <===>  arima.sim(....)

AR <- c(0.7, -0.1, 0.2)
MA <- c(-0.5, 0.4, 0.4)
n <- 512 ; sd <- 0.1
n.st <- 10

set.seed(1)
for(i in 1:200) {
    cat(sprintf("%3d ", i))
    p <- sample(0:length(AR), 1)
    q <- sample(0:length(MA), 1)
    .ar <- AR[seq_len(p)]
    .ma <- MA[seq_len(q)]
    n.st <- p+q+ rpois(1, lambda = 2)
    sid <- round(runif(1)* 1000)
    set.seed(sid)
    y1 <- arima.sim(list(ar = .ar, ma = .ma), n = n, n.start = n.st, sd = sd)
    set.seed(sid)
    y2 <- fracdiff.sim(n = n, ar = .ar, ma = - .ma, d = 0,
                       n.start = n.st, sd = sd, backComp = FALSE)
    if(!isTRUE(aeq <- all.equal(c(y1), y2$series, tol = 1e-15)))
        cat("y1 and y2 are not equal: ", aeq,"\n")
    if(!(i %% 10))    cat("\n")
}


## Last Line:
cat('Time elapsed: ', proc.time(),'\n')


