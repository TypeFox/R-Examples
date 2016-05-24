.jevalIterator <-
function(x)
{
    r <- NULL
    while(.jcall(x, "Z", "hasNext"))
        r <- c(r, list(.jcall(x, "Ljava/lang/Object;", "next")))
    r
}
