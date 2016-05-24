deviance.hdlm <-
function(object, ...) {
    sum(residuals(object)^2, na.rm=TRUE)
}

