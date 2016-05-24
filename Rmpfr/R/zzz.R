## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_Rmpfr_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(sprintf("C code of R package 'Rmpfr': GMP using %d bits per limb\n",
				  .mpfr.gmp.numbbits()))
}

.onLoad <- function(libname, pkgname) {
    if(mpfrVersion() < "3.0.0")
	warning("MPFR C library version ", format(mpfrVersion()),
		" is outdated, and minor functionality will be missing.\n",
		"  Consider installing a newer version of MPFR (e.g., from mpfr.org),\n",
		"  and re-install the R package Rmpfr after that.", call.=FALSE)
}
