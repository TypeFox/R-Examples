setClass(Class = "desirability", representation = representation(response = "character", low = "numeric", high = "numeric", target = "ANY", scale = "numeric", 
    importance = "numeric"))
desirability = function(response, low, high, target = "max", scale = c(1, 1), importance = 1, constraints) {
    if (low >= high) 
        stop("the lower bound must be greater than the high bound!")
    if (any(scale <= 0)) 
        stop("the scale parameter must be greater than zero!")
    if (!is.numeric(target) & !identical(tolower(target), "min") & !identical(tolower(target), "max")) 
        stop("target needs to be \"min\", \"max\" or a numeric value")
    return(new("desirability", response = deparse(substitute(response)), low = low, high = high, target = target, scale = scale, importance = importance))
}
.desireFun = function(low, high, target = "max", scale = c(1, 1), importance = 1) {
    DB = FALSE
    if (importance > 10 | importance < 0.1) 
        stop("importance needs to be in [0.1, 10]")
    if (low >= high) 
        stop("the lower bound must be greater than the high bound!")
    if (any(scale <= 0)) 
        stop("the scale parameter must be greater than zero!")
    if (is.numeric(target)) {
        out = function(y) {
            if (DB) 
                print("target")
            flush.console()
            d = rep(0, length(y))
            d[y >= low & y <= target] = ((y[y >= low & y <= target] - low)/(target - low))^scale[1]
            d[y >= target & y <= high] = ((y[y >= target & y <= high] - high)/(target - high))^scale[2]
            return(d^importance)
        }
        return(out)
    }
    if (identical(tolower(target), "min")) {
        out = function(y) {
            if (DB) 
                print("min")
            d = rep(0, length(y))
            d[y > high] = 0
            d[y < low] = 1
            d[y >= low & y <= high] = ((y[y >= low & y <= high] - high)/(low - high))^scale[1]
            return(d^importance)
        }
        return(out)
    }
    if (identical(tolower(target), "max")) {
        out = function(y) {
            if (DB) 
                print("max")
            d = rep(0, length(y))
            d[y < low] = 0
            d[y > high] = 1
            d[y >= low & y <= high] = ((y[y >= low & y <= high] - low)/(high - low))^scale[1]
            return(d^importance)
        }
        return(out)
    }
}
setMethod("show", signature(object = "desirability"), function(object) {
    if (!is.numeric(object@target)) 
        cat("Target is to", paste(object@target, "imize", sep = ""), object@response, "\n")
    else cat("Target is ", object@target, " for", object@response, "\n")
    cat("lower Bound: ", object@low, "\n")
    cat("higher Bound: ", object@high, "\n")
    if (is.numeric(object@target)) 
        cat("Scale factor is: low =", object@scale[1], "and high =", object@scale[2], "\n")
    else if (identical("min", object@target) | identical("max", object@target)) 
        cat("Scale factor is: ", object@scale, "\n")
    cat("importance: ", object@importance, "\n")
    cat("\n")
})
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "desirability"), function(x, y, scale, main, xlab, ylab, type, col, numPoints = 500, ...) {
    xm1 = NULL
    xm2 = NULL
    ym = NULL
    y = NULL
    if (missing(main)) 
        main = paste("Desirability function for", x@response)
    if (missing(xlab)) 
        xlab = x@response
    if (missing(ylab)) 
        ylab = "Desirability"
    if (missing(type)) 
        type = "l"
    if (missing(scale)) 
        scale = x@scale
    if (missing(col)) 
        col = 1:length(scale)
    dFun = .desireFun(x@low, x@high, x@target, x@scale, x@importance)
    xVals = seq(x@low - 0.04 * diff(range(x@low, x@high)), x@high + 0.04 * diff(range(x@low, x@high)), length = numPoints)
    yVals = dFun(xVals)
    plot(xVals, yVals, main = main, xlab = xlab, ylab = ylab, type = type, col = col, ...)
    if (is.numeric(x@target)) {
        xm1 = mean(c(par("usr")[1], x@target))
        xm2 = mean(c(par("usr")[2], x@target))
        ym1 = yVals[max((1:numPoints)[xVals <= xm1])]
        ym2 = yVals[max((1:numPoints)[xVals <= xm2])]
        text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1, paste("scale =", scale[1]), adj = c(0, 0))
        text(xm2 - 0.025 * diff(range(par("usr")[1:2])), ym2, paste("scale =", scale[2]), adj = c(1, 1))
    }
    else {
        xm1 = mean(par("usr")[c(1, 2)])
        ym1 = yVals[max((1:numPoints)[xVals <= xm1])]
        if (identical(x@target, "max")) 
            text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1 - 0.025 * diff(range(par("usr")[3:4])), paste("scale =", scale[1]), adj = c(0, 0))
        else text(xm1 + 0.025 * diff(range(par("usr")[1:2])), ym1 + 0.025 * diff(range(par("usr")[3:4])), paste("scale =", scale[1]), adj = c(0, 1))
    }
    out = list(x = xVals, y = yVals)
    names(out) = c(x@response, "desirability")
    invisible(out)
})
overall = function(fdo, steps = 20, constraints, ...) {
    DB = FALSE
    importances = list()
    cs = list()
    if (!missing(constraints)) 
        cs = constraints
    l = vector(mode = "list", length = 0)
    fitList = fits(fdo)
    if (length(fitList) < 1) 
        stop(paste("no fits found in fits(", deparse(substitute(fdo)), ")"), sep = "")
    desList = desires(fdo)
    if (length(desList) < 1) 
        stop(paste("no desirabilities found in desires(", deparse(substitute(fdo)), ")"), sep = "")
    X = cube(fdo)
    newdata = NULL
    for (i in names(names(fdo))) {
        seqList = vector(mode = "list", length = 0)
        seqList[["length"]] = steps
        seqList[["from"]] = min(X[, i])
        seqList[["to"]] = max(X[, i])
        minC = NULL
        maxC = NULL
        if (!is.null(cs[[i]])) {
            if (length(cs[[i]]) < 2) 
                stop("length of ", names(cs[i]), "=", cs[i], " < 2 in constraints")
            minC = min(cs[[i]])
            if (!is.null(minC) & !is.na(minC)) 
                seqList[["from"]] = minC
            maxC = max(cs[[i]])
            if (!is.null(maxC) & !is.na(maxC)) 
                seqList[["to"]] = maxC
            if (maxC == minC) 
                stop(paste("equal values in constraints ", names(cs[i]), "=", cs[i]))
        }
        l[[i]] = do.call(seq, seqList)
    }
    if (DB) 
        print(l)
    newdata = expand.grid(l)
    names(newdata) = names(X)
    out = newdata
    yCharSet = intersect(names(desires(fdo)), names(fits(fdo)))
    dFrame = data.frame(matrix(NA, nrow = nrow(newdata), ncol = length(yCharSet) + 1))
    names(dFrame) = c(yCharSet, "overall")
    dFrame[, "overall"] = 1
    for (y in yCharSet) {
        obj = desList[[y]]
        dFun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
        lm.y = fitList[[y]]
        importances[[y]] = desires(fdo)[[y]]@importance
        yHat = predict(lm.y, newdata = newdata, ...)
        yDes = dFun(yHat)
        dFrame[, y] = yDes
        if (DB) {
            print(y)
            print(dFun)
            print(lm.y)
            print(dFrame)
        }
    }
    geomFac = 1/sum(unlist(importances))
    overall = apply(dFrame, 1, prod)^geomFac
    dFrame[, "overall"] = overall
    dFrame = cbind(out, dFrame)
    invisible(dFrame)
}
.desHelp = function(fdo, factors, ...) {
    if (length(factors) != length(names(fdo))) 
        stop("not enough factors specified in factors")
    if (any(is.na(factors))) 
        stop("factors contain NA")
    yCharSet = intersect(names(desires(fdo)), names(fits(fdo)))
    desList = desires(fdo)
    fitList = fits(fdo)
    yDes = vector(mode = "list")
    for (y in yCharSet) {
        obj = desList[[y]]
        dFun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
        lm.y = fitList[[y]]
        yHat = predict(lm.y, newdata = data.frame(factors), ...)
        yDes[[y]] = dFun(yHat)
    }
    return(yDes)
}
setClass(Class = "desOpt", representation = representation(facCoded = "list", facReal = "list", responses = "list", desirabilities = "list", overall = "numeric", 
    all = "data.frame", fdo = "facDesign"))
as.data.frame.desOpt = function(x, row.names = NULL, optional = FALSE, ...) {
    return(x@all)
}
setMethod("as.data.frame", "desOpt", function(x, row.names = NULL, optional = FALSE, ...) {
    return(x@all)
})
.validizeConstraints = function(fdo, constraints) {
    X = as.data.frame(fdo)
    csOut = vector(mode = "list")
    for (i in names(names(fdo))) {
        csOut[[i]] = c(min(X[, i]), max(X[, i]))
    }
    if (missing(constraints)) 
        return(csOut)
    cs2 = constraints[names(names(fdo))]
    cs2 = cs2[!unlist(lapply(cs2, is.null))]
    cs2 = cs2[(unlist(lapply(cs2, length)) == 2)]
    csOut[names(cs2)] = cs2[names(cs2)]
    return(csOut)
}
.dHelp = function(model, dFun) {
    lm1 = model
    d1 = dFun
    out = function(newdata) {
        return(d1(predict(lm1, newdata = newdata)))
    }
    return(out)
}
optimum = function(fdo, constraints, steps = 25, type = "grid", start, ...) {
    DB = FALSE
    if (missing(fdo)) 
        stop("missing fdo!")
    X = as.data.frame(fdo)
    numFac = length(names(fdo))
    if (!(type %in% c("grid", "optim", "gosolnp"))) {
        warning(paste("type =", deparse(substitute(type)), "not found --> using type = \"grid\""))
        type = "grid"
    }
    constraints = .validizeConstraints(fdo, constraints)
    if (missing(start)) 
        start = as.numeric(lapply(constraints, mean))
    lower = numeric(length(constraints))
    upper = numeric(length(constraints))
    for (i in seq(along = constraints)) {
        lower[i] = min(constraints[[i]])
        upper[i] = max(constraints[[i]])
    }
    if (DB) {
        print(constraints)
        print(start)
    }
    desOpt = new("desOpt")
    desOpt@fdo = fdo
    facCoded = NA
    desirabilities = NA
    overall = NA
    setList = list()
    dList = list()
    importances = list()
    yCharSet = intersect(names(desires(fdo)), names(fits(fdo)))
    for (y in yCharSet) {
        obj = desires(fdo)[[y]]
        dFun = .desireFun(obj@low, obj@high, obj@target, obj@scale, obj@importance)
        lm.y = fits(fdo)[[y]]
        importances[[y]] = desires(fdo)[[y]]@importance
        dList[[y]] = .dHelp(lm.y, dFun)
    }
    geomFac = 1/sum(unlist(importances))
    dAll = function(X) {
        newdata = data.frame(t(X))
        names(newdata) = LETTERS[1:ncol(newdata)]
        return(prod(unlist(lapply(dList, do.call, list(newdata = newdata))))^geomFac)
    }
    dAllRsolnp = function(X) {
        newdata = data.frame(t(X))
        names(newdata) = LETTERS[1:ncol(newdata)]
        return(-prod(unlist(lapply(dList, do.call, list(newdata = newdata))))^geomFac)
    }
    if (type == "optim") {
#        print(lower)
 #       print(upper)
        temp = optim(par = start, dAll, method = "L-BFGS-B", lower = lower, upper = upper, control = list(fnscale = -1, maxit = 1000))
        facCoded = as.list(temp$par)
        names(facCoded) = names(names(fdo))
        desOpt@facCoded = facCoded
        overall = temp$value
        desirabilities = .desHelp(fdo, desOpt@facCoded)
    }
    if (type == "gosolnp") {
        #if (!require(Rsolnp, quietly = TRUE)) 
        #    stop("Package Rsolnp needs to be installed!")
        temp = Rsolnp::gosolnp(fun = dAllRsolnp, LB = lower, UB = upper)
        facCoded = as.list(temp$pars)
        names(facCoded) = names(names(fdo))
        desOpt@facCoded = facCoded
        overall = -rev(temp$values)[1]
        desirabilities = .desHelp(fdo, desOpt@facCoded)
    }
    if (type == "grid") {
        dVals = overall(fdo = fdo, constraints = constraints, steps = steps)
        index = order(dVals[, "overall"], decreasing = TRUE)[1]
        desOpt@all = dVals
        desOpt@facCoded = as.list(dVals[index, names(names(fdo))])
        desirabilities = as.list(dVals[index, names(response(fdo))])
        names(desirabilities) = names(response(fdo)) #fix for the case of having just one response 
        overall = dVals[index, "overall"]
    }
    for (i in names(desOpt@facCoded)) {
        desOpt@facReal[[i]] = code2real(lows(fdo)[[i]], highs(fdo)[[i]], desOpt@facCoded[[i]])
    }
    desOpt@desirabilities = desirabilities
    desOpt@overall = overall
    newData = do.call(data.frame, desOpt@facCoded)
    for (i in names(desOpt@desirabilities)) {
        desOpt@responses[[i]] = predict(fits(fdo)[[i]], newData)
    }
    return(desOpt)
}
setMethod("show", signature(object = "desOpt"), function(object) {
    cat(paste("\ncomposite (overall) desirability:", format(object@overall, digits = 3)))
    cat("\n")
    cat("\n")
    temp1 = do.call(data.frame, object@facCoded)
    temp2 = do.call(data.frame, object@facReal)
    facFrame = rbind(temp1, temp2)
    row.names(facFrame) = c("coded", "real")
    show(format(facFrame, digits = 3))
    temp1 = do.call(data.frame, object@responses)
    temp2 = do.call(data.frame, object@desirabilities)
    respDesFrame = rbind(temp1, temp2)
    row.names(respDesFrame) = c("Responses", "Desirabilities")
    cat("\n")
    show(format(respDesFrame, digits = 3))
}) 
