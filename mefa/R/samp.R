samp <- function (x, ...) 
    if (isS4(x)) mefa4::samp(x) else UseMethod("samp")
