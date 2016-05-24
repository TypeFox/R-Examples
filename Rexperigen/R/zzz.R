## We're starting up global state
.onLoad <- function(libname, pkgname){
    op <- options()
    op.setup <- list(
        Rexperigen.experimenter = "",
        Rexperigen.password = "",
        Rexperigen.server = "db.phonologist.org",
        Rexperigen.server.version = "1.0.0"
    )

    ## Thank you Hadley http://r-pkgs.had.co.nz/r.html
    toset <- !(names(op.setup) %in% names(op))
    if(any(toset)) options(op.setup[toset])

  invisible()
}
