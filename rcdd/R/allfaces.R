allfaces <- function(hrep) {

    stopifnot(is.character(hrep) || is.numeric(hrep))
    validcdd(hrep, representation = "H")

    if (is.character(hrep)) {
        .Call("allfaces", hrep, PACKAGE = "rcdd")
    } else {
        storage.mode(hrep) <- "double"
        .Call("allfaces_f", hrep, PACKAGE = "rcdd")
    }
}
