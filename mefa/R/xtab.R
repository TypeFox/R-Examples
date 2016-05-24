xtab <- function (x, ...) 
    if (isS4(x)) mefa4::xtab(x) else UseMethod("xtab")
