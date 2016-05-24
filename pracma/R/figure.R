##
##  f i g u r e . r  Figure (Matlab Style)
##


figure <- function(figno, title = "") {
    if (missing(figno)) figno <- NULL
    else {
        if (!is.numeric(figno) || length(figno) != 1 || floor(figno) != ceiling(figno))
            stop("The figure handle must be a whole number.")
    }

    if (!is.null(figno)) {
        if (figno == 1 || figno == -1) {
            stop("Device 1 is the null device; cannot be opened or closed.")
        } else if (figno == 0) {
            stop("Requested figure handle 0 is in use by another object.")
        }
    }

#     if (.Platform$OS.type == "unix") {
#         if (.Platform$GUI == "AQUA") win <- quartz
#         else if (.Platform$GUI == "X11") win <- X11
#         else
#             stop("Unknown platform GUI for Unix platform type.")
#     } else if (.Platform$OS.type == "Linux") {
#         win <- X11
#     } else if (.Platform$OS.type == "windows") {
#         win <- windows
#     } else
#         stop(paste("Unknown platform type:", .Platform$OS.type))

    if (is.null(figno)) {
        dev.new()               # dev.new() may be platform independent
    } else {
        devl <- dev.list()
        if (abs(figno) %in% devl) {
            if (figno > 0) dev.set(figno)
            else           dev.off(-figno)
        } else {
            cat("Device", figno, "is not an open device. List of devices:\n")
            return(devl)
        }
    }
    invisible()
}
