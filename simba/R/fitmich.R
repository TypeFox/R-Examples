"fitmich" <- 
function(x, y, a=3, b=0.5) {
        res <- nls(y~(a*x)/(b+x), start=list(a=a, b=b))
        res
    }