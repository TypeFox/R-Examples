print.stateCD <-
function(x, ...)
    print(data.frame(a.p=x$a.p, a.s=x$a.s, kc=x$kc, w=x$w))
