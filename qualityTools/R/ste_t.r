setClass(Class = "steepAscent", representation = representation(name = "character", X = "data.frame", response = "data.frame"))
setMethod("response", "steepAscent", function(object) {
    out = object@response
    return(out)
})
setReplaceMethod("response", "steepAscent", function(object, value) {
    if (is.vector(value)) {
        temp = data.frame(value)
        names(temp) = deparse(substitute(value))
        if (nrow(object@X) == nrow(temp)) {
            object@response = temp
            return(object)
        }
        stop("number of rows differ!")
    }
    if (is.data.frame(value)) {
        if (nrow(object@X) == nrow(value)) {
            object@response = value
            return(object)
        }
        stop("number of rows differ!")
    }
    stop(paste(deparse(substitute(value)), " needs to be a vector or data.frame"))
})
setMethod("[", signature(x = "steepAscent", i = "ANY", j = "ANY"), function(x, i, j) {
    bound = ncol(x@X)
    if (j <= bound) 
        x@X[i, j]
    else x@response[i, j - bound]
})
setMethod("as.data.frame", "steepAscent", function(x, row.names = NULL, optional = FALSE, ...) {
    return(cbind(x@X, x@response))
})
as.data.frame.steepAscent = function(x, row.names = NULL, optional = FALSE, ...) {
    return(cbind(x@X, x@response))
}
setMethod("show", signature(object = "steepAscent"), function(object) {
    print(as.data.frame(object))
})
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "steepAscent"), function(x, y, ...) {
    Delta = (x@X)$Delta
    frame = cbind(Delta, response(x))
    names(frame) = c("Delta", names(response(x)))
    plot(frame, ...)
})
steepAscent = function(factors, response, size = 0.2, steps = 5, data) {
    DB = FALSE
    if (missing(data)) 
        return("missing an object of class 'facDesign'")
    else fdo = data
    if (missing(factors) | length(factors) < 1) 
        return("missing factors")
    if (!is.character(factors)) 
        return("factors needs to be a character")
    names(names(fdo))
    model = fits(data)[[response]]
    if (DB) 
        print(model)
    if (is.null(model)) {
        form = c(response, "~")
        for (i in seq(along = factors)) {
            if (i == 1) 
                form = c(form, factors[i])
            else form = c(form, "+", factors[i])
        }
        form = paste(form, collapse = "")
        model = lm(form, data = fdo)
    }
    if (DB) 
        print(model)
    b = numeric(length = length(factors))
    x = numeric(length = length(factors))
    names(x) = factors
    for (i in seq(along = factors)) {
        b[i] = coef(model)[factors[i]]
        if (i == 1) {
            x[i] = size * sign(b[i])
        }
        else {
            if (DB) {
                print(x[1])
                print(b[1])
                print(b[i])
            }
            x[i] = (x[1]/b[1]) * b[i]
        }
    }
    if (DB) 
        print(x)
    Run = 1:(steps + 1)
    Delta = 0:steps
    frameOut = data.frame(Run, Delta)
    initial = ncol(frameOut)
    for (i in seq(along = factors)) {
        frameOut[, i + initial] = x[i] * 0:steps
        names(frameOut)[i + initial] = paste(factors[i], ".coded", collapse = "", sep = "")
        if (DB) 
            print(factors[i])
    }
    initial = ncol(frameOut)
    for (i in seq(along = factors)) {
        frameOut[, i + initial] = code2real(lows(fdo)[[factors[i]]], highs(fdo)[[factors[i]]], x[i] * 0:steps)
        names(frameOut)[i + initial] = paste(factors[i], ".real", collapse = "", sep = "")
        if (DB) 
            print(factors[i])
    }
    soa = new("steepAscent")
    soa@X = frameOut
    response(soa) = rep(NA, times = nrow(frameOut))
    names(response(soa)) = deparse(substitute(response))
    cat("\n")
    cat(paste("Steepest Ascent for", deparse(substitute(data)), "\n"))
    cat("\n")
    print(format(frameOut, digits = 3))
    invisible(soa)
} 
