variable.names.manyglm <- 
function (object, full = FALSE, ...) 
{
    if (full) 
        dimnames(object$qr[[1]]$qr)[[2]]
    else if (object$rank) 
        dimnames(object$qr[[1]]$qr)[[2]][seq_len(object$rank)]
    else character(0)
}

