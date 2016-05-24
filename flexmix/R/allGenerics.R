#
#  Copyright (C) 2004-2012 Friedrich Leisch and Bettina Gruen
#  $Id: allGenerics.R 4834 2012-08-02 10:17:09Z gruen $
#

setGeneric("flexmix",
           function(formula, data=list(), k=NULL,
                    cluster=NULL, model=NULL, concomitant=NULL, control=NULL,
                    weights = NULL)
           standardGeneric("flexmix"))

setGeneric("FLXfit",
           function(model, concomitant, control,
                    postunscaled=NULL, groups, weights)
           standardGeneric("FLXfit"))


###**********************************************************

setGeneric("FLXgetModelmatrix",
           function(model, data, ...) standardGeneric("FLXgetModelmatrix"))

setGeneric("FLXfillConcomitant",
           function(concomitant, ...) standardGeneric("FLXfillConcomitant"))

###**********************************************************

setGeneric("logLik")

setGeneric("clogLik", function(object, ...) standardGeneric("clogLik"))

setGeneric("EIC", function(object, ...) standardGeneric("EIC"))

###**********************************************************

setGeneric("FLXcheckComponent", function(model, ...) standardGeneric("FLXcheckComponent"))

setGeneric("FLXgetK", function(model, ...) standardGeneric("FLXgetK"))

setGeneric("FLXgetObs", function(model) standardGeneric("FLXgetObs"))

setGeneric("FLXmstep", function(model, ...) standardGeneric("FLXmstep"))

setGeneric("FLXremoveComponent", function(model, ...) standardGeneric("FLXremoveComponent"))

setGeneric("FLXdeterminePostunscaled", function(model, ...) standardGeneric("FLXdeterminePostunscaled"))

setGeneric("FLXgetDesign", function(object, ...) standardGeneric("FLXgetDesign"))

setGeneric("FLXreplaceParameters", function(object, ...) standardGeneric("FLXreplaceParameters"))

setGeneric("FLXlogLikfun", function(object, ...) standardGeneric("FLXlogLikfun"))

setGeneric("FLXgradlogLikfun", function(object, ...) standardGeneric("FLXgradlogLikfun"))

setGeneric("VarianceCovariance", function(object, ...) standardGeneric("VarianceCovariance"))

setGeneric("FLXgetParameters", function(object, ...) standardGeneric("FLXgetParameters"))

setGeneric("logLikfun_comp", function(object, ...) standardGeneric("logLikfun_comp"))

setGeneric("getPriors", function(object, ...) standardGeneric("getPriors"))

setGeneric("existGradient", function(object, ...) standardGeneric("existGradient"))

setGeneric("refit_mstep", function(object, newdata, ...)  standardGeneric("refit_mstep"))
setGeneric("refit_optim", function(object, ...)  standardGeneric("refit_optim"))

###**********************************************************

setGeneric("group", function(object, ...) standardGeneric("group"))
setGeneric("rflexmix", function(object, newdata, ...) standardGeneric("rflexmix"))
setGeneric("rFLXM", function(model, components, ...) standardGeneric("rFLXM"))

## just to make sure that some S3 generics are available in S4
setGeneric("fitted")
setGeneric("predict")
setGeneric("simulate")
setGeneric("summary")
setGeneric("unique")
