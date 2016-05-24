# real77: A typeII + scaling exponent (assuming replacement)
# Now fully deprecated

real77 <- function(...) {
    msg <- paste0("The real77* family was deprecated in version 0.5.\n",
                  "A more sensible alternative (the flexp* family) has been added.\n",
                  "See ?real77 for more information.\n")
    stop(msg, call. = FALSE)  
}

real77_fit <- function(...) {
    real77(...)
}	

real77_nll <- function(...) {
    real77(...)
}	

real77_diff <- function(...) {
    real77(...)
}	

real77_nll_diff <- function(...) {
    real77(...)
}	

