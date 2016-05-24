print.trk <- function(x, ...) {
    cat("\t", "Stellar track\n\n")
    cat("Mass =", x$mass,"Msun\n")
    cat("Z =", x$z, ", Y =", x$y, "\n")
    cat("Mixing length =", x$ml, "\n")
    cat("[alpha/Fe] =", x$alpha.enh, "\n")
    invisible(x)
}

print.zahb <- function(x, ...) {
    cat("\t", "ZAHB\n\n")
    RGB <- x$data$mass[length(x$data$mass)]
    cat("RGB progenitor mass =", round(RGB, 2), "Msun \n")
    cat("Z =", x$z, ", Y =", x$y, "\n")
    cat("Mixing length =", x$ml, "\n")
    cat("[alpha/Fe] =", x$alpha.enh, "\n")
    invisible(x)
}

print.iso <- function(x, ...) {
    cat("\t", "Stellar isochrone\n\n")
    cat("Age =", x$age,"Gyr\n")
    cat("Z =", x$z, ", Y =", x$y, "\n")
    cat("Mixing length =", x$ml, "\n")
    cat("[alpha/Fe] =", x$alpha.enh, "\n")
    invisible(x)
}

print.hb <- function(x, ...) {
    cat("\t", "Stellar track from ZAHB\n\n")
    cat("Mass =", x$mass,"Msun\n")
    cat("Mass RGB =", x$massRGB,"Msun\n")
    cat("Z =", x$z, ", Y =", x$y, "\n")
    cat("Mixing length =", x$ml, "\n")
    cat("[alpha/Fe] =", x$alpha.enh, "\n")
    invisible(x)
}
