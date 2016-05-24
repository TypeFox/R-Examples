
# $Id: AAA.R 270 2006-09-07 11:28:44Z hothorn $

.onLoad <- function(lib, pkg) {
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
    return(TRUE)
}
