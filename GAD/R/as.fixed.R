as.fixed <-
function (x)
{
f <- factor(x)
class(f) <- c("factor", "fixed")
f
}

