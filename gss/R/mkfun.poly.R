## Make RK for cubic splines
mkrk.cubic <- function(range)
{
    ## Create the environment
    env <- list(min=min(range), max=max(range))
    ## Create the rk function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ##% Check the inputs
        if (!(is.vector(x)&is.vector(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        if ((min(x,y)<env$min)|(max(x,y)>env$max)) {
            stop("gss error in rk: inputs are out of range")
        }
        ##% Scale the inputs
        x <- (x-env$min)/(env$max-env$min)
        y <- (y-env$min)/(env$max-env$min)
        ##% Return the result
        rk <- function(x,y) {
            k2 <- function(x) ((x-.5)^2-1/12)/2
            k4 <- function(x) ((x-.5)^4-(x-.5)^2/2+7/240)/24
            k2(x)*k2(y)-k4(abs(x-y))
        }
        if (outer.prod) outer(x,y,rk)
        else rk(x,y)
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make phi function for cubic splines
mkphi.cubic <- function(range)
{
    ## Create the environment
    env <- list(min=min(range), max=max(range))
    ## Create the phi function
    fun <- function(x,nu,env) {
        ##% Check the input
        if (!is.vector(x)) {
            stop("gss error in phi: inputs are of wrong types")
        }
        if ((min(x)<env$min)|(max(x)>env$max)) {
            stop("gss error in phi: inputs are out of range")
        }
        ##% Return the result
        (x-env$min)/(env$max-env$min)-.5
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make RK for periodic cubic splines
mkrk.cubic.per <- function(range)
{
    ## Create the environment
    env <- list(min=min(range), max=max(range))
    ## Create the rk function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ##% Check the inputs
        if (!(is.vector(x)&is.vector(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        if ((min(x,y)<env$min)|(max(x,y)>env$max)) {
            stop("gss error in rk: inputs are out of range")
        }
        ##% Scale the inputs
        x <- (x-env$min)/(env$max-env$min)
        y <- (y-env$min)/(env$max-env$min)
        ##% Return the result
        rk <- function(x,y) {
            k4 <- function(x) ((x-.5)^4-(x-.5)^2/2+7/240)/24
            -k4(abs(x-y))
        }
        if (outer.prod) outer(x,y,rk)
        else rk(x,y)
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make RK for trigonometric splines
mkrk.trig <- function(range)
{
    ## Create the environment
    env <- list(min=min(range), max=max(range))
    ## Create the rk function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ##% Check the inputs
        if (!(is.vector(x)&is.vector(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        if ((min(x,y)<env$min)|(max(x,y)>env$max)) {
            stop("gss error in rk: inputs are out of range")
        }
        ##% Scale the inputs
        x <- (x-env$min)/(env$max-env$min)
        y <- (y-env$min)/(env$max-env$min)
        ##% Return the result
        rk <- function(x,y) {
            k4 <- function(x) ((x-.5)^4-(x-.5)^2/2+7/240)/24
            -k4(abs(x-y))-2*cos(2*pi*(x-y))/(2*pi)^4
        }
        if (outer.prod) outer(x,y,rk)
        else rk(x,y)
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make phi function for trigonometric splines
mkphi.trig <- function(range)
{
    ## Create the environment
    env <- list(min=min(range), max=max(range))
    ## Create the phi function
    fun <- function(x,nu,env) {
        ##% Check the input
        if (!is.vector(x)) {
            stop("gss error in phi: inputs are of wrong types")
        }
        if ((min(x)<env$min)|(max(x)>env$max)) {
            stop("gss error in phi: inputs are out of range")
        }
        ##% Return the result
        xx <- (x-env$min)/(env$max-env$min)
        switch(nu,cos(2*pi*xx),sin(2*pi*xx))
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make RK for linear splines
mkrk.linear <- function(range)
{
    ## Create the environment
    env <- list(min=min(range), max=max(range))
    ## Create the rk function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ##% Check the inputs
        if (!(is.vector(x)&is.vector(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        if ((min(x,y)<env$min)|(max(x,y)>env$max)) {
            stop("gss error in rk: inputs are out of range")
        }
        ##% Scale the inputs
        x <- (x-env$min)/(env$max-env$min)
        y <- (y-env$min)/(env$max-env$min)
        ##% Return the result
        rk <- function(x,y) {
            k1 <- function(x) (x-.5)
            k2 <- function(x) ((x-.5)^2-1/12)/2
            k1(x)*k1(y)+k2(abs(x-y))
        }
        if (outer.prod) outer(x,y,rk)
        else rk(x,y)
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make RK for periodic linear splines
mkrk.linear.per <- function(range)
{
    ## Create the environment
    env <- list(min=min(range), max=max(range))
    ## Create the rk function
    fun <- function(x,y,env,outer.prod=FALSE) {
        ##% Check the inputs
        if (!(is.vector(x)&is.vector(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        if ((min(x,y)<env$min)|(max(x,y)>env$max)) {
            stop("gss error in rk: inputs are out of range")
        }
        ##% Scale the inputs
        x <- (x-env$min)/(env$max-env$min)
        y <- (y-env$min)/(env$max-env$min)
        ##% Return the result
        rk <- function(x,y) {
            k2 <- function(x) ((x-.5)^2-1/12)/2
            k2(abs(x-y))
        }
        if (outer.prod) outer(x,y,rk)
        else rk(x,y)
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}
