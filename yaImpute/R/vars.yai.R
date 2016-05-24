# vars returns a list of variables used in the problem
# xvars and yvars only returns x's and y's
xvars = function (object) colnames(object$xRefs)
yvars = function (object) colnames(object$yRefs)
vars  = function (object) list(xvars=colnames(object$xRefs),yvars=colnames(object$yRefs))
