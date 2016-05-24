##
##  b e e p . R
##


beep <- function() cat("\a")


disp <- function(...) cat(..., "\n")


ver <- function() {
    z <- list()
    z$R.version <- R.Version()
    z$platform <- z$R.version$platform
    if (nzchar(.Platform$r_arch)) 
        z$platform <- paste(z$platform, .Platform$r_arch, sep = "/")
    z$platform <- paste(z$platform, " (", 8 * .Machine$sizeof.pointer, 
        "-bit)", sep = "")
    z$locale <- Sys.getlocale()

    cat("------------------------------------------------------------------------\n")
    cat("Version:    ", z$R.version$version.string, "\n", sep = "")
    cat("License:    Gnu General Public License, GPL-3\n", sep = "")
    cat("Platform:   ", z$R.version$platform, "\n", sep = "")
    cat("Op.System:  ",  z$platform, "\n", sep = "")
    cat("Locale:     ",  z$locale, "\n", sep = "")
    cat("------------------------------------------------------------------------\n")

    # Loaded base and contributed packages
    session_info <- sessionInfo()
    base_pkgs <- session_info$basePkgs
    cat("Loaded Base Packages:\n  ", base_pkgs, "\n", sep = "  ")
    cat("Contributed Packages:\n")
    other_pkgs <- session_info$otherPkgs
    for (pack in other_pkgs) {
        pkg_name <- pack$Package; l <- max(12 - nchar(pkg_name), 2)
        pkg_name <- paste(pkg_name, blanks(l), sep = "", collapse = "")

        pkg_version <- pack$Version; l <- max(8 - nchar(pkg_version), 2)
        pkg_version <- paste("Version ", pkg_version, blanks(l), sep = "", collapse = "")

        pkg_date <- substr(pack$Date, 1, 10)
        pkg_date <- paste("(", pkg_date, ")", sep = "", collapse = "")

        pkg_license <- pack$License
        pkg_license <- paste("  License ", pkg_license, sep = "", collapse = "")

        cat(pkg_name, pkg_version, pkg_date, pkg_license, "\n", sep = "  ")
    }
    cat("------------------------------------------------------------------------\n")
    invisible(NULL)
}
