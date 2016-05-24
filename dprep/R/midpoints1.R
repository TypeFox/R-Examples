midpoints1 <-
function (x) {
    n <- length(x)
mid=rep(0,n)
points <- .C("Points", as.double(x), as.integer(n), mpoint = as.double(mid),PACKAGE="dprep")
points$mpoint
}
