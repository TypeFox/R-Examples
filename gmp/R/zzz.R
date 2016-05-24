## we "need" S4 methods for dispatch on both (x,y)  .noGenerics <- TRUE
.conflicts.OK <- TRUE

.gmpVersion <- function() .Call(R_gmp_get_version)
gmpVersion <- function()
    numeric_version(sub("^([0-9]+\\.[0-9]+\\.[0-9]+).*","\\1", .gmpVersion()))

.onLoad <- function(libname, pkgname) {
    options("gmp:warnModMismatch" = TRUE, ## see ../man/biginteger.Rd
            "gmp:warnNoInv" = TRUE) ## ../man/add.biginteger.Rd | ../src/bigmod.cc
}

