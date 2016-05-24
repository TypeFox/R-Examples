setClass("MSALinearity", representation(X = "data.frame", Y = "data.frame", model = "lm", conf.level = "numeric", Linearity = "numeric", GageName = "character", 
    GageTolerance = "numeric", DateOfStudy = "character", PersonResponsible = "character", Comments = "character", facNames = "character"))
setMethod("response", "MSALinearity", function(object) {
    out = object@Y
    return(out)
})
setReplaceMethod("response", "MSALinearity", function(object, value) {
    if (is.vector(value) == TRUE) 
        value = data.frame(matrix(value, ncol = ncol(object@Y)))
    if (is.matrix(value) == TRUE) 
        value = data.frame(value)
    object@Y = value
    return(object)
})
setMethod("summary", signature(object = "MSALinearity"), function(object) {
    cat("----------------------", fill = TRUE)
    print(object)
    cat("----------------------", fill = TRUE)
    print(summary(object@model))
    cat("----------------------", fill = TRUE)
    names(object@Linearity) = "Linearity:"
    print(object@Linearity)
})
setMethod("plot", signature(x = "MSALinearity"), function(x, ylim, col, pch, lty = c(1, 2), ...) {
    if (is(x, "MSALinearity") == FALSE) 
        stop("object needs to be from class 'MSALinearity'")
    conf.level = x@conf.level
    g = nrow(x@X[2])
    m = ncol(x@Y)
    if (missing(col)) 
        col = c(1, 2, 1, 4)
    if (missing(pch)) 
        pch = c(20, 18)
    bias = x@Y
    mbias = numeric(g)
    for (i in 1:g) bias[i, ] = x@Y[i, ] - x@X[i, 2]
    for (i in 1:g) mbias[i] = mean(as.numeric(bias[i, ]))
    if (missing(ylim)) 
        ylim = c(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE))
    plot(x = x@X$Ref, y = mbias, ylim = ylim, col = col[2], pch = pch[2], ylab = "Bias", xlab = "Reference Values", ...)
    for (i in 1:g) points(x = rep(x@X$Ref[i], length = m), y = bias[i, ], col = col[1], pch = pch[1])
    abline(h = 0, lty = 3, col = "gray")
    BIAS = numeric()
    ref = numeric()
    for (i in 1:g) {
        BIAS = c(BIAS, as.numeric(bias[i, ]))
        ref = c(ref, rep(x@X$Ref[i], length = m))
    }
    lm.1 = lm(formula = BIAS ~ ref)
    a = lm.1[[1]][2]
    names(a) = "slope"
    b = lm.1[[1]][1]
    names(b) = "intercept"
    y.vec = numeric()
    for (i in 1:g) y.vec = c(y.vec, x@Y[i, ])
    pre = predict.lm(lm.1, interval = "confidence", level = conf.level)
    lines(ref, pre[, 1], col = col[3], lty = lty[1])
    lines(ref, pre[, 2], col = col[4], lty = lty[2])
    lines(ref, pre[, 3], col = col[4], lty = lty[2])
    legend("topright", legend = c("Single Bias", "Mean Bias", "Regression", paste(conf.level * 100, "% conf.level")), pch = c(pch, -1, -1), col = col, lty = c(-1, 
        -1, lty), inset = 0.04)
})
setMethod("show", signature(object = "MSALinearity"), function(object) {
    print(as.data.frame(object))
})
as.data.frame.MSALinearity = function(x, row.names = NULL, optional = FALSE, ...) {
    return(cbind(x@X, x@Y))
}
setMethod("as.data.frame", "MSALinearity", function(x, row.names = NULL, optional = FALSE, ...) {
    return(cbind(x@X, x@Y))
})
gageLin <- function(object, conf.level = 0.95, ylim, col, pch, lty = c(1, 2), stats = TRUE, plot = TRUE) {
    if (is(object, "MSALinearity") == FALSE) 
        stop("object needs to be from class 'MSALinearity'")
    object@conf.level = conf.level
    g = nrow(object@X[2])
    m = ncol(object@Y)
    if (missing(col)) 
        col = c(1, 2, 1, 4)
    if (missing(pch)) 
        pch = c(20, 18)
    bias = object@Y
    mbias = numeric(g)
    for (i in 1:g) bias[i, ] = object@Y[i, ] - object@X[i, 2]
    for (i in 1:g) mbias[i] = mean(as.numeric(bias[i, ]))
    if (missing(ylim)) 
        ylim = c(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE))
    if (plot == TRUE) 
        plot(x = object@X$Ref, y = mbias, ylim = ylim, col = col[2], pch = pch[2], ylab = "Bias", xlab = "Reference Values")
    if (plot == TRUE) {
        for (i in 1:g) points(x = rep(object@X$Ref[i], length = m), y = bias[i, ], col = col[1], pch = pch[1])
        abline(h = 0, lty = 3, col = "gray")
    }
    BIAS = numeric()
    ref = numeric()
    for (i in 1:g) {
        BIAS = c(BIAS, as.numeric(bias[i, ]))
        ref = c(ref, rep(object@X$Ref[i], length = m))
    }
    lm.1 = lm(formula = BIAS ~ ref)
    object@model <- lm.1
    a = lm.1[[1]][2]
    names(a) = "slope"
    b = lm.1[[1]][1]
    names(b) = "intercept"
    y.vec = numeric()
    for (i in 1:g) y.vec = c(y.vec, object@Y[i, ])
    if (plot == TRUE) {
        pre = predict.lm(lm.1, interval = "confidence", level = conf.level)
        lines(ref, pre[, 1], col = col[3], lty = lty[1])
        lines(ref, pre[, 2], col = col[4], lty = lty[2])
        lines(ref, pre[, 3], col = col[4], lty = lty[2])
    }
    test = numeric(g + 1)
    for (i in 1:g) test[i] = t.test(bias[, i], mu = 0, conf.level = conf.level)["p.value"]
    test[g + 1] = t.test(BIAS, mu = 0, conf.level = conf.level)["p.value"]
    Linearity = abs(lm.1[[1]][2]) * 100
    object@Linearity = Linearity
    names(Linearity) = "LINEARITY:"
    if (plot == TRUE) {
        legend("topright", legend = c("Single Bias", "Mean Bias", "Regression", paste(conf.level * 100, "% conf.level")), pch = c(pch, -1, -1), col = col, lty = c(-1, 
            -1, lty), inset = 0.04)
    }
    if (stats == TRUE) {
        cat("----------------------", fill = TRUE)
        cat("BIAS:", fill = TRUE)
        print(bias)
        cat("----------------------", fill = TRUE)
        cat("MEAN OF BIAS:", fill = TRUE)
        temp = mbias
        names(temp) = rownames(object@Y)
        print(temp)
        cat("----------------------", fill = TRUE)
        cat("LINEAR MODEL:", fill = TRUE)
        print(summary(lm.1))
        cat("----------------------", fill = TRUE)
        print(Linearity)
        return(object)
    }
    else return(object)
}
setClass("MSALinearity", representation(X = "data.frame", Y = "data.frame", model = "lm", conf.level = "numeric", Linearity = "numeric", GageName = "character", 
    GageTolerance = "numeric", DateOfStudy = "character", PersonResponsible = "character", Comments = "character", facNames = "character"))
setMethod("response", "MSALinearity", function(object) {
    out = object@Y
    return(out)
})
setReplaceMethod("response", "MSALinearity", function(object, value) {
    if (is.vector(value) == TRUE) 
        value = data.frame(matrix(value, ncol = ncol(object@Y)))
    if (is.matrix(value) == TRUE) 
        value = data.frame(value)
    object@Y = value
    return(object)
})
setMethod("summary", signature(object = "MSALinearity"), function(object) {
    cat("----------------------", fill = TRUE)
    print(object)
    cat("----------------------", fill = TRUE)
    print(summary(object@model))
    cat("----------------------", fill = TRUE)
    names(object@Linearity) = "Linearity:"
    print(object@Linearity)
})
setMethod("plot", signature(x = "MSALinearity"), function(x, ylim, col, pch, lty = c(1, 2), ...) {
    if (is(x, "MSALinearity") == FALSE) 
        stop("object needs to be from class 'MSALinearity'")
    conf.level = x@conf.level
    g = nrow(x@X[2])
    m = ncol(x@Y)
    if (missing(col)) 
        col = c(1, 2, 1, 4)
    if (missing(pch)) 
        pch = c(20, 18)
    bias = x@Y
    mbias = numeric(g)
    for (i in 1:g) bias[i, ] = x@Y[i, ] - x@X[i, 2]
    for (i in 1:g) mbias[i] = mean(as.numeric(bias[i, ]))
    if (missing(ylim)) 
        ylim = c(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE))
    plot(x = x@X$Ref, y = mbias, ylim = ylim, col = col[2], pch = pch[2], ylab = "Bias", xlab = "Reference Values", ...)
    for (i in 1:g) points(x = rep(x@X$Ref[i], length = m), y = bias[i, ], col = col[1], pch = pch[1])
    abline(h = 0, lty = 3, col = "gray")
    BIAS = numeric()
    ref = numeric()
    for (i in 1:g) {
        BIAS = c(BIAS, as.numeric(bias[i, ]))
        ref = c(ref, rep(x@X$Ref[i], length = m))
    }
    lm.1 = lm(formula = BIAS ~ ref)
    a = lm.1[[1]][2]
    names(a) = "slope"
    b = lm.1[[1]][1]
    names(b) = "intercept"
    y.vec = numeric()
    for (i in 1:g) y.vec = c(y.vec, x@Y[i, ])
    pre = predict.lm(lm.1, interval = "confidence", level = conf.level)
    lines(ref, pre[, 1], col = col[3], lty = lty[1])
    lines(ref, pre[, 2], col = col[4], lty = lty[2])
    lines(ref, pre[, 3], col = col[4], lty = lty[2])
    legend("topright", legend = c("Single Bias", "Mean Bias", "Regression", paste(conf.level * 100, "% conf.level")), pch = c(pch, -1, -1), col = col, lty = c(-1, 
        -1, lty), inset = 0.04)
})
setMethod("show", signature(object = "MSALinearity"), function(object) {
    print(as.data.frame(object))
})
as.data.frame.MSALinearity = function(x, row.names = NULL, optional = FALSE, ...) {
    return(cbind(x@X, x@Y))
}
setMethod("as.data.frame", "MSALinearity", function(x, row.names = NULL, optional = FALSE, ...) {
    return(cbind(x@X, x@Y))
})
gageLinDesign <- function(ref, n = 5) {
    numColY = n
    numRowY = length(ref)
    X = data.frame(cbind(Part = 1:length(ref), Ref = ref))
    Y = data.frame(matrix(NA, nrow = numRowY, ncol = numColY))
    for (i in 1:n) names(Y)[i] = paste("n", i, sep = "")
    gageLinDesign = new("MSALinearity")
    gageLinDesign@X = X
    gageLinDesign@Y = Y
    return(gageLinDesign)
} 
