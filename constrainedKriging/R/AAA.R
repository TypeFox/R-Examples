.CONSTRAINEDKRIGING_CACHE <- new.env(FALSE, parent=globalenv())

.onAttach <- function(lib, pkg) {
    assign("gpclib", FALSE, envir=.CONSTRAINEDKRIGING_CACHE)
    packageStartupMessage(paste("\nThe constrainedKriging package provides functions for efficient",
                                "\ncomputations of nonlinear spatial predictions with local change of support.\n\n"
	      ))
}

.onUnload <- function(libpath) {
    rm(.CONSTRAINEDKRIGING_CACHE)
}

