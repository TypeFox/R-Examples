
"update.CPTgrain" <- function(object,  ...){

    if(!(object$isCompiled))
        object <- compile( object )

    ##cl <- match.call(expand.dots=TRUE)
    args <- list(...)
    arg.names <- names(args)

    if ("CPTlist" %in% arg.names){
        object$cptlist[names(args$CPTlist)] <- args$CPTlist
        pot.with.1        <- .createPotList( object$rip, object$universe )
        newpot            <- .insertCPT(object$cptlist, pot.with.1, details=0)
        object$origpot    <- newpot
        object$temppot    <- newpot
        object$equipot    <- .insertNA(pot.with.1)
        object$isPropagated <- FALSE
    }
    object
}





















