`MC` <-
function(f, env = NULL) if(is.R()) f else
{
    if(is.R()) return(f)
    if (mode(f) != "function") stop(paste("not a function:", f))
    if (length(env) > 0 && any(names(env <- as.list(env)) == ""))
        stop(paste("contains unnamed arguments:", env))
    fargs <- if (length(f) > 1) f[1:(length(f) - 1)] else NULL
    if (any(duplicated(names(fargs <- c(fargs,env)))))
        stop(paste("duplicated arguments:", paste(names(fargs)),
            collapse = ", "))
    cf <- as.function(c(fargs, f[length(f)]))
    return(cf) }
