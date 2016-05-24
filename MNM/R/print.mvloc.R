`print.mvloc` <-
function(x,...)
    {
    X<-list(location=x$location,vcov=x$vcov)
    print(X, ...)
    }
