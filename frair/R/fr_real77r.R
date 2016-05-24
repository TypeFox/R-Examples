# real77r: A typeII + scaling exponent (not assuming replacement)
# Now fully deprecated

real77r <- function(...) {
    msg <- paste0("The real77* family was deprecated in version 0.5.\n",
                  "A more sensible alternative (the flexp* family) has been added.\n",
                  "See ?real77 for more information.\n")
    stop(msg, call. = FALSE)  
}

real77r_fit <- function(...) {
    real77(...)
}	

real77r_nll <- function(...) {
    real77(...)
}	

real77r_diff <- function(...) {
    real77(...)
}	

real77r_nll_diff <- function(...) {
    real77(...)
}	

