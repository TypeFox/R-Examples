##### Chargement de base



## An environment useful to store the options
.adehabitatMAEnv <- new.env()

.onLoad <- function(lib, pkg)
{
    environment(.adehabitatMAEnv) <- asNamespace("adehabitatMA")
    assign(".adeoptions", list(epsilon=1e-8),
           envir=.adehabitatMAEnv)

}




adeoptions <- function(...)
{
    olde <- get(".adeoptions", envir=.adehabitatMAEnv)
    class(olde) <- "optade"
    oo <- list(...)
    if (length(oo)>0) {
        if (is.list(oo[[1]])) {
            if (inherits(oo[[1]], "optade"))
                oo <- oo[[1]]
        }
        newe <- olde
        for (i in names(oo)) {
            newe[[i]] <- oo[[i]]
        }
        assign(".adeoptions", newe, envir=.adehabitatMAEnv)
        invisible(olde)
    } else {
        return(olde)
    }
}

