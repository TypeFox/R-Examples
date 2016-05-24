## this prints the summary
setMethod("summary", "pva", function(object) {
    y <- summary(as(object, "dcmle"), fancyPVAmodel(object))
    y@coef <- object@summary
    show(y)
})

## this extracts mcmc info
## and transforms it to original scale if desired
diagn_scale <-
function(object) {
    m <- as(as(object, "dcmle"), "MCMClist")
    object@model@transf(m)
 }

## coef method (takes into account fixed values)
setMethod("coef", "pva", function(object) {
    object@fullcoef
})

## need to define this as well
setMethod("vcov", "pva", function(object) {
    vc <- object@vcov
    fx <- object@model@fixed
    if (is.null(fx))
        return(vc)
    vn <- object@model@varnames
    rv <- matrix(NA, length(vn), length(vn))
    dimnames(rv) <- list(vn, vn)
    i <- !(vn %in% names(fx))
    rv[i,i] <- vc
    rv
})
setMethod("confint", "pva", function(object, diagn_scale=FALSE) {
    if (diagn_scale) {
        objt <- diagn_scale(object)
    } else {
        objt <- as(as(object, "dcmle"), "MCMClist")
    }
    ci <- confint(objt)
    fx <- object@model@fixed
    if (is.null(fx))
        return(ci)
    vn <- object@model@varnames
    vno <- varnames(object)
    vnt <- varnames(objt)
    vnf <- vnt[match(vn, vno)]
    vnf[is.na(vnf)] <- vn[is.na(vnf)]
    rv <- ci[match(vnf, rownames(ci)),]
    rownames(rv)[is.na(rownames(rv))] <- vnf[is.na(rownames(rv))]
    rv
})
