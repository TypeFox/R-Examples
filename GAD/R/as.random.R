as.random <-
function (x) 
{
f <- factor(x)
class(f) <- c("factor", "random")
f
}

