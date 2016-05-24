############### ufn ####################
# function defined in order to deal with out of bounds functions/parameters
ufn <- function(par, fnuser) {
    userfn<-fnuser$fn # Because do.call doesn't know about fnuser
    nlist<-length(fnuser$dots)
    lnames<-names(fnuser$dots)
    dots<-fnuser$dots
    parps<-par*fnuser$PARSCALE
    testf<-try(tryf<-do.call("userfn",c(list(parps),dots)))
    # Note: Idea from R_inferno of using c(list(), dots), but ?? no names(dots) %in% spec
    # Since function may NOT always use par, do not list(par=par).
    if ((class(testf) == "try-error") || is.na(tryf) || is.null(tryf) || 
        is.infinite(tryf)) {
        tryf <- .Machine$double.xmax
        attr(tryf, "inadmissible") <- TRUE
    }
    else {
        attr(tryf, "inadmissible") <- FALSE
    }
    if (is.null(tryf)) stop("NULL FUNCTION")
    fnuser$KFN<-1+fnuser$KFN
    if (fnuser$MAXIMIZE) tryf <- -tryf # handle the maximization
    tryf/fnuser$FNSCALE # and scale to finish
}
############## end ufn ###################


