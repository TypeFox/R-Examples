ci.pois.pearson.hartley.approx <-
function (x, alpha = 0.05, ci.type = c("two-sided", "lower", 
    "upper")) 
{
    n <- length(x)
    x <- sum(x)
    ci.type <- match.arg(ci.type)
    if (x == 0) 
        stop("All values of 'x' are 0.  You must use the exact method in this case.")
    vl <- 2 * x
    vu <- 2 * (x + 1)
    switch(ci.type, `two-sided` = {
        lcl <- qchisq(alpha/2, vl)/2
        ucl <- qchisq(1 - (alpha/2), vu)/2
    }, lower = {
        lcl <- qchisq(alpha, vl)/2
        ucl <- Inf
    }, upper = {
        lcl <- 0
        ucl <- qchisq(1 - alpha, vu)/2
    })
    ci.limits <- c(lcl, ucl)/n
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = "lambda", 
        limits = ci.limits, type = ci.type, method = "Pearson/Hartley approx", 
        conf.level = 1 - alpha, sample.size = n)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
