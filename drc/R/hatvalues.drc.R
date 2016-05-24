hatvalues.drc <- function(model, ...)
{
    xmat <- model$der
    diag(xmat %*% ginv(t(xmat) %*% xmat) %*% t(xmat))
#    names(hvector) <- as.character(1:length(hvector))
#    hvector
}

