### To find the boundary of a bounded convex set
"convex.bounds" <-
function (x, dir, indFunc, ..., tol = 1e-07) 
{
    ## x: a point within the set
    ## dir: a vector giving the direction along which bounds are sought
    ## indFunc: the indicator function of a bounded convex set
    ## ... : additional arguments passed to indFunc
    if (all(dir == 0)) 
        stop("invalid direction in convex.bounds()")
    if (indFunc(x, ...) < 0.5) 
        stop("x not in the support of indFunc")
    f.onedim <- function(u) indFunc(x + u * dir, ...)
    e <- -2
    while (f.onedim(e) > 0.5) e <- e * 2
    lower <- e
    e <- 2
    while (f.onedim(e) > 0.5) e <- e * 2
    upper <- e
    ans <- numeric(2)
    ## search for `lower' boundary along dir
    bracket.low <- lower
    bracket.high <- 0
    repeat {
        cand <- 0.5 * (bracket.low + bracket.high)
        if (f.onedim(cand) > 0.5) 
            bracket.high <- cand
        else bracket.low <- cand
        if (bracket.high - bracket.low < tol) {
            ans[1] <- bracket.high
            break
        }
    }
    ## search for `upper' boundary along dir
    bracket.low <- 0
    bracket.high <- upper
    repeat {
        cand <- 0.5 * (bracket.low + bracket.high)
        if (f.onedim(cand) > 0.5) 
            bracket.low <- cand
        else bracket.high <- cand
        if (bracket.high - bracket.low < tol) {
            ans[2] <- bracket.low
            break
        }
    }
    return(ans)
}

### Wrapper to arms.c
"arms" <-
function (y.start, myldens, indFunc, n.sample, ...) 
{
    ## y.start: starting point
    ## myldens: univariate or multivariate logdensity from which a sample
    ##          needs to be generated
    ## indFunc: the indicator function of the support of myldens
    ##          (assumed to be convex and bounded)
    ## n.sample: desired sample size
    ## ...     : additional arguments passed to myldens and indFunc
    ## sanity checks first
#     if (mode(myldens) != "function") 
#         stop("myldens not a function")
#     if (mode(indFunc) != "function") 
#         stop("indFunc not a function")
#     if (n.sample < 0) 
#         stop("n.sample must be nonnegative")
#     if (n.sample < 1) 
#         return(numeric(0))
#     if (!is.numeric(y.start)) 
#         stop("non numeric argument y.start")
    dim <- length(y.start)
#     if (dim == 0) 
#         stop("starting point has length zero")
#     if (!(indFunc(y.start, ...) > 0)) 
#         stop("starting point not in the support")
    if (dim == 1) {
        bounds <- y.start + convex.bounds(y.start, dir = 1, indFunc, 
            ...)
        if ( diff(bounds) < 1e-7 )
            y.sample <- rep(y.start, n.sample)
        else {
            f <- function(x) myldens(x, ...)
            y.sample <- .Call("arms", bounds, f, y.start, as.integer(n.sample), 
                              new.env())
        }
    }
    else {
        y.sample <- rbind(y.start, matrix(0, n.sample, dim))
        for (k in 1:n.sample) {
            ## pick a direction at random 
            dir <- rnorm(dim)
            ## look for boundaries of support in the selected direction
            bounds <- convex.bounds(y.sample[k, ], dir, indFunc, 
                ...)
            if ( diff(bounds) < 1e-7 )
                y.sample[k + 1, ] <- y.sample[k, ]
            else {
                ## define the univariate density to be passed to arms.c
                f <- function(x) myldens(y.sample[k, ] + x * dir, 
                                         ...)
                ## call arms.c
                y.sample[k + 1, ] <- y.sample[k, ] + dir * .Call("arms", 
                    bounds, f, 0, as.integer(1), new.env())
            }
        }
        y.sample <- y.sample[-1, ]
    }
    return(y.sample)
}

