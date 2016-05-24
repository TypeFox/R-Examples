taxa <- function (x, ...) 
    if (isS4(x)) mefa4::taxa(x) else UseMethod("taxa")
