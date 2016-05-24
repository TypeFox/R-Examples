svTaskCallbackManager <- function (handlers = list(), registered = FALSE,
verbose = FALSE) 
{
    suspended <- FALSE
    .verbose <- verbose
    
    add <- function(f, data = NULL, name = NULL, register = TRUE) {
        if (is.null(name)) 
            name <- as.character(length(handlers) + 1L)
        handlers[[name]] <<- list(f = f)
        if (!missing(data)) 
            handlers[[name]][["data"]] <<- data
        if (!registered && register) {
            register()
        }
        name
    }
    
    remove <- function(which) {
        if (is.character(which)) {
            tmp <- (1L:length(handlers))[!is.na(match(which, 
                names(handlers)))]
            if (length(tmp)) 
                stop(gettextf("no such element '%s'", which), 
                  domain = NA)
            which <- tmp
        }
        else which <- as.integer(which)
        handlers <<- handlers[-which]
        return(TRUE)
    }
    
    evaluate <- function(expr, value, ok, visible) {
        if (suspended) 
            return(TRUE)
        discard <- character(0L)
        for (i in names(handlers)) {
            h <- handlers[[i]]
            if (length(h) > 1L) {
                val <- h[["f"]](expr, value, ok, visible, i[["data"]])
            }
            else {
                val <- h[["f"]](expr, value, ok, visible)
            }
            if (!val) {
                discard <- c(discard, i)
            }
        }
        if (length(discard)) {
            if (.verbose) 
                cat(gettext("Removing"), paste(discard, collapse = ", "), 
                  "\n")
            idx <- is.na(match(names(handlers), discard))
            if (length(idx)) 
                handlers <<- handlers[idx]
            else handlers <<- list()
        }
        return(TRUE)
    }
    
    suspend <- function(status = TRUE) {
        suspended <<- status
    }
    
    register <- function(name = "SV-taskCallbackManager", verbose = .verbose) {
        if (verbose) 
            cat(gettext("Registering evaluate as low-level callback\n"))
        id <- addTaskCallback(evaluate, name = name)
        registered <<- TRUE
        id
    }
    
    list(add = add, evaluate = evaluate, remove = remove, register = register, 
        suspend = suspend, callbacks = function() handlers)
}
