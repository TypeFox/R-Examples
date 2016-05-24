## The "incorrect" implimentation of the Beddington-DeAngelis Type-II curve
# Now fully deprecated 

bdII <- function(...) {
    msg <- paste0("The bdII model was deprecated in version 0.5.\n",
                  "Use bdaII or emdII depending on your needs.\n",
                  "See ?bdII for more information.\n")
    stop(msg, call. = FALSE)  
}

bdII_fit <- function(...) {
    bdII(...)
}	

bdII_nll <- function(...) {
    bdII(...) 
}

bdII_diff <- function(...) {
    bdII(...) 
}

bdII_nll_diff <- function(...) {
    bdII(...)
}