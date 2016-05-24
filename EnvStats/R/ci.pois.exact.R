ci.pois.exact <-
function (x, alpha = 0.05, ci.type = c("two-sided", "lower", 
    "upper")) 
{
    n <- length(x)
    x <- sum(x)
    ci.type <- match.arg(ci.type)
    if (x == 0) {
        if (ci.type != "upper") {
            ci.type <- "upper"
            warning("All values of 'x' are 0; 'ci.type' forced to be 'upper'")
        }
        lcl <- 0
        ucl <- -log(alpha)
    }
    else {
        fcn.to.optimize.lcl <- function(lambda, x, p) {
            ((1 - ppois(x - 1, lambda)) - p)^2
        }
        fcn.to.optimize.ucl <- function(lambda, x, p) {
            (ppois(x, lambda) - p)^2
        }
        ci.obj <- ci.pois.pearson.hartley.approx(x, alpha = alpha, 
            ci.type = ci.type)$limits
        switch(ci.type, `two-sided` = {
            lcl <- nlminb(start = ci.obj["LCL"], objective = fcn.to.optimize.lcl, 
                lower = .Machine$double.eps, x = x, p = alpha/2)$par
            ucl <- nlminb(start = ci.obj["UCL"], objective = fcn.to.optimize.ucl, 
                lower = .Machine$double.eps, x = x, p = alpha/2)$par
        }, lower = {
            lcl <- nlminb(start = ci.obj["LCL"], objective = fcn.to.optimize.lcl, 
                lower = .Machine$double.eps, x = x, p = alpha)$par
            ucl <- Inf
        }, upper = {
            lcl <- 0
            ucl <- nlminb(start = ci.obj["UCL"], objective = fcn.to.optimize.ucl, 
                lower = .Machine$double.eps, x = x, p = alpha)$par
        })
    }
    ci.limits <- c(lcl, ucl)/n
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = "lambda", 
        limits = ci.limits, type = ci.type, method = "exact", 
        conf.level = 1 - alpha, sample.size = n)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
