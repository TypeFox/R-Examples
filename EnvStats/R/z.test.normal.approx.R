z.test.normal.approx <-
function (theta.hat, sd.theta.hat, hyp.theta, alternative) 
{
    z.stat <- (theta.hat - hyp.theta)/sd.theta.hat
    names(z.stat) <- "z"
    p.value <- switch(alternative, two.sided = {
        2 * (1 - pnorm(abs(z.stat)))
    }, greater = {
        1 - pnorm(z.stat)
    }, less = {
        pnorm(z.stat)
    })
    ret.obj <- list(statistic = z.stat, p.value = p.value, estimate = c(theta.hat = theta.hat), 
        null.value = c(theta = hyp.theta), alternative = alternative, 
        method = "One-sample z-test")
    ret.obj
}
