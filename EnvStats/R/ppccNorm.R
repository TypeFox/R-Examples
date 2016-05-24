ppccNorm <-
function (x) 
{
    x <- sort(x)
    m.tilda <- qnorm(ppoints(x, a = 3/8))
    ss.m <- sum(m.tilda^2)
    c.vec <- m.tilda/sqrt(ss.m)
    cor(c.vec, x)
}
