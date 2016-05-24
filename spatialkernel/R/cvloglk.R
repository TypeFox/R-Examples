cvloglk <- function(pts, marks, t = NULL, h)
{
  if(is.null(t)) {
    ans <- cvlogl(pts, marks, h)
  }	else ans <- cvloglp(pts, marks, t, h)
  ans$hcv <- h[which.max(ans$cv)]
  ans
}

## pts[n*2], y[n].
## h[] can be unequally separated values
cvlogl <- function(pts, marks, h)
{
    ##if(exists(".adaptpara", env=.GlobalEnv)) {
    ##    .adaptpara <- get(".adaptpara", env=.GlobalEnv) ##, env=.GlobalEnv)
    ##} else .adaptpara <- get(".adaptpara", env = getNamespace("spatialkernel"))
	adapt <- chkernel()
    n <- length(marks)
    nh <- length(h)
    types <- unique(marks)
    mtypes <- 1:length(types) - 1 ## y must from 0 to m-1
    names(mtypes) <- types
    y <- mtypes[marks]
    c <- NULL
    for(i in 1:nh) c <- cbind(c, rep(1, n))
    ans<-.C("lcn", as.double(pts), as.integer(y), as.integer(n), as.double(h), 
        as.integer(nh), as.integer(adapt$kernel), as.double(c),
        lc=double(nh), PACKAGE="spatialkernel")$lc
    invisible(list(cv=ans, pts=pts, marks=marks, h=h))
}

## pooled cvlogl
cvloglp <- function(pts, marks, t, h)
{
    tt <- sort(unique(t))
    ntt <- length(tt)
    lcp <- rep(0, length(h))
    for(i in 1:ntt) {
        ndx <- which(t==tt[i])
        lcp <- lcp+cvlogl(pts[ndx,], marks[ndx], h)$cv
    }
    invisible(list(cv=lcp, pts=pts, marks=marks, t=t, h=h))
}
