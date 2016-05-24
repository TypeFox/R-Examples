## Convert a byte-compiled function to an interpreted-code function
unByteCode <- function(fun)
{
    FUN <- eval(parse(text=deparse(fun)))
    environment(FUN) <- environment(fun)
    FUN
}

## Replace function definition inside of a locked environment **HACK**
assignEdgewise <- function(name, env, value)
{
    unlockBinding(name, env=env)
    assign( name, envir=env, value=value)
    lockBinding(name, env=env)
    invisible(value)
}

## Replace byte-compiled function in a locked environment with an
## interpreted-code function
unByteCodeAssign <- function(fun)
{
    name <- gsub('^.*::+','', deparse(substitute(fun)))
    FUN <- unByteCode(fun)
    retval <- assignEdgewise(name=name,
                             env=environment(FUN),
                             value=FUN
                             )
    invisible(retval)
}
