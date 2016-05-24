ci.binom.exact <-
function (x, size, alpha = 0.05, ci.type = c("two-sided", "lower", 
    "upper")) 
{
    if (length(x) != 1 || length(size) != 1) 
        stop("x and size must be vectors of length 1")
    ci.type <- match.arg(ci.type)
    if (x == 0 || x == size) {
        if (x == 0) {
            if (ci.type != "upper") {
                ci.type <- "upper"
                warning("All values of 'x' are 0; 'ci.type' forced to be 'upper'")
            }
            lcl <- 0
            ucl <- 1 - alpha^(1/size)
        }
        else {
            if (ci.type != "lower") {
                ci.type <- "lower"
                warning("All values of 'x' are 1; 'ci.type' forced to be 'lower'")
            }
            lcl <- alpha^(1/size)
            ucl <- 1
        }
    }
    else {
        p.hat <- x/size
        v1 <- 2 * (size - x + 1)
        v2 <- 2 * x
        switch(ci.type, `two-sided` = {
            lcl <- x/(x + (size - x + 1) * qf(1 - alpha/2, v1, 
                v2))
            num <- (x + 1) * qf(1 - alpha/2, v2 + 2, v1 - 2)
            ucl <- num/(size - x + num)
        }, lower = {
            lcl <- x/(x + (size - x + 1) * qf(1 - alpha, v1, 
                v2))
            ucl <- 1
        }, upper = {
            lcl <- 0
            num <- (x + 1) * qf(1 - alpha, v2 + 2, v1 - 2)
            ucl <- num/(size - x + num)
        })
    }
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = "prob", 
        limits = ci.limits, type = ci.type, method = "Exact (Clopper-Pearson)", 
        conf.level = 1 - alpha, sample.size = size)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
