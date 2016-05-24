#
#  Copyright (C) 2004-2013 Friedrich Leisch and Bettina Gruen
#  $Id: concomitant.R 4922 2013-09-03 13:32:45Z gruen $
#

FLXPmultinom <- function(formula=~1) {
  z <- new("FLXPmultinom", name="FLXPmultinom", formula=formula)
  multinom.fit <- function(x, y, w, ...) {
    r <- ncol(x)
    p <- ncol(y)
    if (p < 2) stop("Multinom requires at least two components.")
    mask <- c(rep(0, r + 1), rep(c(0, rep(1, r)), p - 1))
    nnet::nnet.default(x, y, w, mask = mask, size = 0, 
                       skip = TRUE, softmax = TRUE, censored = FALSE, 
                       rang = 0, trace=FALSE,...)
  }
  z@fit <- function(x, y, w, ...) multinom.fit(x, y, w, ...)$fitted.values
  z@refit <- function(x, y, w, ...) {
    if (missing(w) || is.null(w)) w <- rep(1, nrow(y))
    fit <- multinom(y ~ 0 + x, weights = w, data = list(y = y, x = x), Hess = TRUE, trace = FALSE)
    fit$coefnames <- colnames(x)
    fit$vcoefnames <- fit$coefnames[seq_along(fit$coefnames)]
    dimnames(fit$Hessian) <- lapply(dim(fit$Hessian) / ncol(x), function(i) paste(rep(seq_len(i) + 1, each = ncol(x)), colnames(x), sep = ":"))
    fit
  }
  z
}

FLXPconstant <- function() {
  new("FLXP", name="FLXPconstant", formula = ~1,
      fit = function(x, y, w, ...){
        if (missing(w) || is.null(w)) return(matrix(colMeans(y), ncol=ncol(y), dimnames = list("prior", seq_len(ncol(y)))))
        else return(matrix(colMeans(w*y)/mean(w), ncol=ncol(y), dimnames = list("prior", seq_len(ncol(y)))))
      })
}

###**********************************************************

setMethod("FLXgetModelmatrix", signature(model="FLXP"),
function(model, data, groups, lhs, ...)
{
    mt <- terms(model@formula, data=data)
    mf <- model.frame(delete.response(mt), data=data, na.action = NULL)
    X <- model.matrix(mt, data=mf)
    if (nrow(X)){
        if (!checkGroup(X, groups$group))
            stop("model variables have to be constant for grouping variable")
        model@x <- X[groups$groupfirst,,drop=FALSE]
    }
    else{
        model@x <- matrix(1, nrow=sum(groups$groupfirst))
    }
  
    model
})

checkGroup <- function(x, group) {
  check <- TRUE
  for(g in levels(group)){
    gok <- group==g
    if(any(gok)){
      check <- all(c(check, apply(x[gok,,drop=FALSE], 2, function(z) length(unique(z)) == 1)))
    }
  }
  check
}

###**********************************************************

setMethod("refit_mstep", signature(object="FLXP", newdata="missing"),
function(object, newdata, posterior, group, ...) NULL)

setMethod("refit_mstep", signature(object="FLXPmultinom", newdata="missing"),
function(object, newdata, posterior, group, ...) {
  groupfirst <- if (length(group)) groupFirst(group) else rep(TRUE, nrow(posterior))
  object@refit(object@x, posterior[groupfirst,,drop=FALSE], ...)
})

###**********************************************************


setMethod("FLXfillConcomitant", signature(concomitant="FLXP"), function(concomitant, posterior, weights) {
  concomitant@coef <- concomitant@refit(concomitant@x, posterior, weights)
  concomitant
})

setMethod("FLXfillConcomitant", signature(concomitant="FLXPmultinom"), function(concomitant, posterior, weights) {
  concomitant@coef <- cbind("1" = 0, t(coef(concomitant@refit(concomitant@x, posterior, weights))))
  concomitant
})


###**********************************************************

