tTestNormal.approx <-
function (theta.hat, sd.theta.hat, hyp.theta, df, alternative) 
{
    t.stat <- (theta.hat - hyp.theta)/sd.theta.hat
    names(t.stat) <- "t"
    p.value <- switch(alternative, two.sided = {
        2 * (1 - pt(abs(t.stat), df = df))
    }, greater = {
        1 - pt(t.stat, df = df)
    }, less = {
        pt(t.stat, df = df)
    })
    ret.obj <- list(statistic = t.stat, parameters = c(df = df), 
        p.value = p.value, estimate = c(theta.hat = theta.hat), 
        null.value = c(theta = hyp.theta), alternative = alternative, 
        method = "One-sample t-Test")
    ret.obj
}
