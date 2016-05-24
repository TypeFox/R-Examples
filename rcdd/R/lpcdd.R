lpcdd <- function(hrep, objgrd, objcon = as(0, class(objgrd)), minimize = TRUE,
    solver = c("DualSimplex", "CrissCross")) {

    solver <- match.arg(solver)
    stopifnot(is.character(hrep) || is.numeric(hrep))
    stopifnot(is.character(objgrd) || is.numeric(objgrd))
    stopifnot(is.character(objcon) || is.numeric(objcon))
    stopifnot(is.character(hrep) == is.character(objgrd))
    stopifnot(is.character(hrep) == is.character(objcon))
    stopifnot(is.logical(minimize))

    stopifnot(ncol(hrep) - 2 == length(objgrd))
    stopifnot(length(objcon) == 1)
    stopifnot(length(minimize) == 1)

    validcdd(hrep, representation = "H")

    if (is.character(hrep)) {
        .Call("lpcdd", hrep, c(objcon, objgrd), minimize, solver,
            PACKAGE = "rcdd")
    } else {
        storage.mode(hrep) <- "double"
        objgrd <- as.double(objgrd)
        objcon <- as.double(objcon)
        .Call("lpcdd_f", hrep, c(objcon, objgrd), minimize, solver,
            PACKAGE = "rcdd")
    }
}
