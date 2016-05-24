.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Please use predefined Credentials only for the testing requests. To obtain your own Credentials see help(authorize).")
}

.onLoad <- function(libname, pkgname) {
    op <- options()
    op.rga <- list(
        rga.client.id = "144394141628-8m5i5icva7akegi3tp6215d9eg9o5cln.apps.googleusercontent.com",
        rga.client.secret = "wlFmhluHqTdZw6UG22h5A2nr",
        rga.cache = ".ga-token.rds"
    )
    toset <- !(names(op.rga) %in% names(op))
    if (any(toset)) options(op.rga[toset])

    invisible()
}
