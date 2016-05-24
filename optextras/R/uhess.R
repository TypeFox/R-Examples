############### uhess.R ####################
# ?? need tests of scaling to make sure we have everything right
uhess <- function(par, fnuser) {
    npar <- length(par)
    uhess<-fnuser$hess
    parps<-par*fnuser$PARSCALE
    dots<-fnuser$dots
    th <- try(tryh <- do.call("uhess", c(list(parps), dots)), silent = TRUE)
    if ((class(th) == "try-error") || any(is.na(tryh)) || any(is.null(tryh)) || 
        any(is.infinite(tryh))) {
        tryh <- matrix(.Machine$double.xmax, nrow = npar, ncol = npar)
        attr(tryh, "inadmissible") <- TRUE
        cat("INADMISSIBLE\n")
    } else {
        attr(tryh, "inadmissible") <- FALSE # assume OK until proven otherwise
    }
    tattr<-attributes(tryh) # save all attrubutes
    if (any(is.null(tryh))) stop("NULL FUNCTION")
    fnuser$KHESS<-1+fnuser$KHESS
    if (fnuser$MAXIMIZE) tryh <- -tryh # handle the maximization
    tryh<- diag(fnuser$PARSCALE)%*%(tryh)%*%diag(fnuser$PARSCALE) 
    # attributes NOT inherited in operation above
    tryh<-tryh/fnuser$FNSCALE
    attributes(tryh)<-tattr # restore attributes (?? not counting ihess??)
    tryh
}  # end uhess definition
############# end uhess ##########################
