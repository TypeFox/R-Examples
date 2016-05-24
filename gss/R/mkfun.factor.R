## Make RK for nominal shrinkage
mkrk.nominal <- function(levels)
{
    k <- length(levels)
    if (k<2) stop("gss error: factor should have at least two levels")
    code <- 1:k
    names(code) <- as.character(levels)
    ## Create the environment
    env <- list(code=code,table=diag(k)-1/k)
    ## Create the rk function
    fun <- function(x, y, env, outer.prod = FALSE) {
        if (!(is.factor(x)&is.factor(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        x <- as.numeric(env$code[as.character(x)])
        y <- as.numeric(env$code[as.character(y)])
        if (any(is.na(c(x,y)))) {
            stop("gss error in rk: unknown factor levels")
        }
        if (outer.prod) env$table[x, y]
        else env$table[cbind(x,y)]
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}

## Make RK for ordinal shrinkage
mkrk.ordinal <- function(levels)
{
    k <- length(levels)
    if (k<2) stop("gss error: factor should have at least two levels")
    code <- 1:k
    names(code) <- as.character(levels)
    ## penalty matrix
    if (k==2) {
        B <- diag(.25,2)
        B[1,2] <- B[2,1] <- -.25
    }
    else {
        B <- diag(2,k)
        B[1,1] <- B[k,k] <- 1
        diag(B[-1,-k]) <- diag(B[-k,-1]) <- -1
        ## Moore-Penrose inverse
        B <- eigen(B)
        B <- B$vec[,-k] %*% diag(1/B$val[-k]) %*% t(B$vec[,-k])
        tol <- sqrt(.Machine$double.eps)
        B <- ifelse(abs(B)<tol,0,B)
    }
    ## Create the environment
    env <- list(code=code,table=B)
    ## Create the rk function
    fun <- function(x, y, env, outer.prod = FALSE) {
        if (!(is.factor(x)&is.factor(y))) {
            stop("gss error in rk: inputs are of wrong types")
        }
        x <- as.numeric(env$code[as.character(x)])
        y <- as.numeric(env$code[as.character(y)])
        if (any(is.na(c(x,y)))) {
            stop("gss error in rk: unknown factor levels")
        }
        if (outer.prod) env$table[x, y]
        else env$table[cbind(x,y)]
    }
    ## Return the function and the environment
    list(fun=fun,env=env)
}
