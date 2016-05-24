# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.







attrassignlm <- function(lmobj) {
	attrassign(model.matrix(lmobj),terms(lmobj))
}

attrassigndefault <- function(mmat,tt) {
        if (!inherits(tt,"terms"))
                stop("need terms object")
        aa<-attr(mmat,"assign")
        if (is.null(aa))
                stop("argument is not really a model matrix")
        ll<-attr(tt,"term.labels")
        if (attr(tt,"intercept")>0)
                ll<-c("(Intercept)",ll)
        aaa<-factor(aa,labels=ll)
        split(order(aa),aaa)
}


if (!isGeneric("attrassign"))
    setGeneric("attrassign", function(object, ...)
        standardGeneric("attrassign"))

setMethod("attrassign", "lm",
         function(object, ...)
         attrassignlm(object, ...))



