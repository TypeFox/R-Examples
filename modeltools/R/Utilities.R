MEapply <- function(object, FUN, clone = TRUE, ...) standardGeneric(MEapply)

setMethod("MEapply", "ModelEnv",
function(object, FUN, clone = TRUE, ...)
{
    ## If we check here, we don't have to check for the existence
    ## of hook collections every time
    if(is.null(FUN))
        return(object)
    
    z <- object
    if (clone)
        z <- clone(object, copydata = FALSE)

    for (name in ls(object@env)){
        if(is.list(FUN)){
            if(name %in% names(FUN)){
                assign(name,
                       FUN[[name]](object@get(name), ...),
                       envir = z@env)
            }
        } else {
            assign(name,
                   FUN(object@get(name), ...),
                   envir = z@env)
        }
    }

    return(z)
})

