## The only guaranteed scores-like item are the projected coordinates
## of 'points'. Usually we get only that, but there may be other
## scores-like things and I anticipate these items may change in the
## future, and we prepare to get any scores-like object.
`scores.ordiplot3d` <-
    function(x, display, ...)
{
    if (missing(display))
        return(x$points)
    ## not the default of points: see what it could be and return
    scoreslike <- names(x)[sapply(x, is.matrix)]
    display <- match.arg(display, c(scoreslike, "sites"))
    ## vegan standards say that scores(x, "sites") should work with
    ## all scores() methods: return points.
    if (display == "sites")
        display <- "points"
    x[[display]]
}
