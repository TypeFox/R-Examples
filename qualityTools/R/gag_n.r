setClass("gageRR", representation = representation(X = "data.frame", ANOVA = "aov", RedANOVA = "aov", method = "character", Estimates = "list", Varcomp = "list", 
    Sigma = "numeric", GageName = "character", GageTolerance = "numeric", DateOfStudy = "character", PersonResponsible = "character", Comments = "character", 
    b = "factor", a = "factor", y = "numeric", facNames = "character", numO = "numeric", numP = "numeric", numM = "numeric"))
setMethod("show", signature(object = "gageRR"), function(object) {
    print(as.data.frame(object))
})
setMethod("[", signature(x = "gageRR", i = "ANY", j = "ANY"), function(x, i, j) {
    x@X[i, j]
})
setMethod("summary", signature(object = "gageRR"), function(object) {
    if (all(is.na(object@X$Measurement))) 
        return(print(as.data.frame(object)))
    cat("\n")
    cat(paste("Operators:\t", object@numO, "\tParts:\t", object@numP))
    cat("\n")
    cat(paste("Measurements:\t", object@numM, "\tTotal:\t", nrow(object@X)))
    cat("\n")
    cat("----------")
    cat("\n")
    return(gageRR(object, method = object@method))
})
setMethod("response", "gageRR", function(object) {
    out = object@X$Measurement
    return(out)
})
setReplaceMethod("response", "gageRR", function(object, value) {
    object@X$Measurement = value
    return(object)
})
setMethod("names", signature(x = "gageRR"), function(x) {
    return(names(as.data.frame(x)))
})
setMethod("as.data.frame", "gageRR", function(x, row.names = NULL, optional = FALSE, ...) {
    return(x@X)
})
as.data.frame.gageRR = function(x, row.names = NULL, optional = FALSE, ...) {
    return(x@X)
}
setGeneric("tolerance", function(x) standardGeneric("tolerance"))
setGeneric("tolerance<-", function(x, value) standardGeneric("tolerance<-"))
setMethod("tolerance", "gageRR", function(x) unlist(x@GageTolerance))
setReplaceMethod("tolerance", "gageRR", function(x, value) {
    if (!is.numeric(value)) 
        stop(paste(deparse(substitute(value)), "needs to be numeric"))
    x@GageTolerance = value
    return(x)
})
setGeneric("sigma", function(x) standardGeneric("sigma"))
setGeneric("sigma<-", function(x, value) standardGeneric("sigma<-"))
setMethod("sigma", "gageRR", function(x) unlist(x@Sigma))
setReplaceMethod("sigma", "gageRR", function(x, value) {
    if (!is.numeric(value)) 
        stop(paste(deparse(substitute(value)), "needs to be numeric"))
    x@Sigma = value
    return(x)
})
.aip = function(x.factor, trace.factor, response, fun = mean, type = c("l", "p", "b"), legend = FALSE, trace.label = deparse(substitute(trace.factor)), 
    fixed = FALSE, xlab = deparse(substitute(x.factor)), ylab = ylabel, ylim = range(cellNew, na.rm = TRUE), lty = nc:1, col = 1, pch = c(1L:9, 0, letters), 
    xpd = NULL, leg.bg = par("bg"), leg.bty = "o", xtick = FALSE, xaxt = par("xaxt"), axes = TRUE, title = "", ...) {
    ylabel <- paste(deparse(substitute(fun)), "of ", deparse(substitute(response)))
    type <- match.arg(type)
    cellNew <- tapply(response, list(x.factor, trace.factor), fun)
    nr <- nrow(cellNew)
    nc <- ncol(cellNew)
    xvals <- 1L:nr
    if (is.ordered(x.factor)) {
        wn <- getOption("warn")
        options(warn = -1)
        xnm <- as.numeric(levels(x.factor))
        options(warn = wn)
        if (!any(is.na(xnm))) 
            xvals <- xnm
    }
    xlabs <- rownames(cellNew)
    ylabs <- colnames(cellNew)
    nch <- max(sapply(ylabs, nchar, type = "width"))
    if (is.null(xlabs)) 
        xlabs <- as.character(xvals)
    if (is.null(ylabs)) 
        ylabs <- as.character(1L:nc)
    xlim <- range(xvals)
    xleg <- xlim[2L] + 0.05 * diff(xlim)
    xlim <- xlim + c(-0.2/nr, if (legend) 0.2 + 0.02 * nch else 0.2/nr) * diff(xlim)
    matplot(xvals, cellNew, ..., type = type, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, xaxt = "n", col = col, lty = lty, pch = pch)
    if (axes && xaxt != "n") {
        axisInt <- function(x, main, sub, lwd, bg, log, asp, ...) axis(1, x, ...)
        mgp. <- par("mgp")
        if (!xtick) 
            mgp.[2L] <- 0
        axisInt(1, at = xvals, labels = xlabs, tick = xtick, mgp = mgp., xaxt = xaxt, ...)
    }
    if (legend) {
        legend("topright", legend = ylabs, title = title, col = col, pch = if (type %in% c("p", "b")) 
            pch, lty = if (type %in% c("l", "b")) 
            lty, bty = leg.bty, bg = leg.bg, inset = 0.02)
    }
    legend("topright", legend = ylabs, title = title, col = col, pch = if (type %in% c("p", "b")) 
        pch, lty = if (type %in% c("l", "b")) 
        lty, bty = leg.bty, bg = leg.bg, inset = c(-0.2, 0), xpd = TRUE)
    invisible()
}
gageRRDesign = function(Operators = 3, Parts = 10, Measurements = 3, method = "crossed", sigma = 6, randomize = TRUE) {
    method = method
    opvec = factor()
    partvec = factor()
    yName = aName = bName = abName = NA
    yName = "Measurement"
    aName = "Operator"
    bName = "Part"
    abName = "Operator:Part"
    Operators = unique(Operators)
    Parts = unique(Parts)
    if (!is.numeric(sigma)) 
        stop("sigma needs to be numeric")
    if (method != "nested" && method != "crossed") 
        warning("unknown method specified --> defaulting to \"method = crossed\"")
    if (!is.numeric(Measurements)) 
        stop("Number of Measurements per Part not specified!")
    else Measurements = round(Measurements[1])
    if (!is.numeric(Operators) && !is.character(Operators)) 
        stop("Operator needs to be numeric 'Operator = 3' or character 'Operator = c(\"A\",\"B\", \"C\"")
    if (is.numeric(Operators)) 
        opvec = factor(LETTERS[1:Operators[1]])
    if (is.character(Operators)) 
        opvec = factor(Operators)
    if (length(unique(opvec)) > 26) 
        stop("To many Operators!")
    if (length(unique(opvec)) < 2) 
        stop("Not enough Operators")
    if (!is.numeric(Parts) && !is.character(Parts)) 
        stop("Parts needs to be numeric 'Parts = 3' or character 'Parts = c(\"A\",\"B\", \"C\"")
    if (is.numeric(Parts)) 
        partvec = factor(LETTERS[1:Parts[1]])
    if (is.character(Parts)) 
        partvec = factor(Parts)
    if (length(unique(partvec)) > 26) 
        stop("To many Parts!")
    if (length(unique(partvec)) < 2) 
        stop("To few Parts")
    Measurement = rep(NA, (length(opvec) * length(partvec) * Measurements))
    outFrame = data.frame()
    if (method == "crossed") {
        temp = expand.grid(opvec, partvec)
        o = rep(temp[, 1], Measurements)
        p = rep(temp[, 2], Measurements)
    }
    else {
        p = rep(sort(rep(partvec, length(opvec))), Measurements)
        o = (rep(opvec, length(Measurement)/length(opvec)))
        p = p[order(o,p)]
        o = o[order(o,p)]
    }
    if (randomize) 
        outFrame = data.frame(StandardOrder = 1:length(Measurement), RunOrder = sample(1:length(Measurement), length(Measurement)), Operator = factor(o), Part = factor(p), 
            Measurement)
    else outFrame = data.frame(StandardOrder = 1:length(Measurement), RunOrder = 1:length(Measurement), Operator = factor(o), Part = factor(p), Measurement)
    outFrame = outFrame[order(outFrame$RunOrder), ]
    gageRRObj = new("gageRR")
    gageRRObj@facNames = c(yName, aName, bName, abName)
    names(gageRRObj@facNames) = c("yName", "aName", "bName", "abName")
    gageRRObj@Sigma = sigma
    gageRRObj@method = method
    gageRRObj@a = factor(o)
    gageRRObj@b = factor(p)
    gageRRObj@y = as.numeric(Measurement)
    gageRRObj@method = method
    gageRRObj@Sigma = sigma
    gageRRObj@X = outFrame
    return(gageRRObj)
}
gageRR = function(gdo, method = "crossed", sigma = 6, alpha = 0.25, DM = NULL, HM = NULL, tolerance = NULL, dig = 3, ...) {
    if (method %in% c("crossed", "nested")) 
        method = method
    else method = gdo@method
    yName = names(gdo)[5]
    aName = names(gdo)[3]
    bName = names(gdo)[4]
    if(method == "crossed")
      abName = paste(aName, ":", bName, sep = "")
    if(method == "nested")
      abName = paste(bName, "(", aName, ")", sep = "")
    bTobName = paste(bName, "to", bName, sep = " ")
    a = gdo@X[, aName]
    b = gdo@X[, bName]
    y = gdo@X[, yName]
    nestedFormula = as.formula(paste(yName, "~", aName, "/", bName))
    crossedFormula = as.formula(paste(yName, "~", aName, "*", bName))
    reducedFormula = as.formula(paste(yName, "~", aName, "+", bName))
    if (!is.null(tolerance)) 
        tolerance(gdo) = tolerance
    if (is.na(y) || !is.numeric(y)) 
        stop("Measurements need to be numeric")
    if (method == "nested") {
        numA <- nlevels(a[, drop = T])
        numB <- nlevels(b[, drop = T])
        numMPP <- length(y)/((numB) * numA)
        gdo@numO = numA
        gdo@numP = numB
        gdo@numM = numMPP
        fit = aov(nestedFormula, data = gdo)
        meanSq <- anova(fit)[, 3]
        gdo@ANOVA = fit
        gdo@method = "nested"
        MSa = meanSq[1]
        MSab = meanSq[2]
        MSe = meanSq[3]
        Cerror = MSe
        Cb = (MSab - MSe)/numMPP
        Ca = (MSa - MSab)/(numB * numMPP)
        if (Ca <= 0) 
            Ca = 0
        if (Cb <= 0) 
            Cb = 0
        Cab = 0
        totalRR = Ca + Cab + Cerror
        repeatability = Cerror
        reproducibility = Ca
        bTob = Cb
        totalVar = Cb + Ca + Cab + Cerror
        estimates = list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
        varcomp = list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, bTob = bTob, totalVar = totalVar)
        gdo@Estimates = estimates
        gdo@Varcomp = varcomp
    }
    if (method == "crossed") {
        numA <- nlevels(a[, drop = T])
        numB <- nlevels(b[, drop = T])
        numMPP <- length(a)/(numA * numB)
        gdo@numO = numA
        gdo@numP = numB
        gdo@numM = numMPP
        fit = aov(crossedFormula, data = gdo)
        model <- anova(fit)
        gdo@ANOVA = fit
        gdo@method = "crossed"
        MSb = MSa = MSab = MSe = 0
        if (bName %in% row.names(model)) 
            MSb = model[bName, "Mean Sq"]
        else warning(paste("missing factor", bName, "in model"))
        if (aName %in% row.names(model)) 
            MSa = model[aName, "Mean Sq"]
        else warning(paste("missing factor", aName, "in model"))
        if (abName %in% row.names(model)) 
            MSab = model[abName, "Mean Sq"]
        else warning(paste("missing interaction", abName, "in model"))
        if ("Residuals" %in% row.names(model)) 
            MSe = model["Residuals", "Mean Sq"]
        else warning("missing Residuals in model")
        Cb = Ca = Cab = Cerror = 0
        Cb = (MSb - MSab)/(numA * numMPP)
        Ca = (MSa - MSab)/(numB * numMPP)
        Cab = (MSab - MSe)/(numMPP)
        Cerror = (MSe)
        gdo@RedANOVA = gdo@ANOVA
        if ((Cab < 0) || (model[abName, "Pr(>F)"] >= alpha)) {
            redFit <- aov(reducedFormula, data = gdo)
            model <- anova(redFit)
            MSb = MSa = MSab = MSe = 0
            if (bName %in% row.names(model)) 
                MSb = model[bName, "Mean Sq"]
            else warning(paste("missing factor", bName, "in model"))
            if (aName %in% row.names(model)) 
                MSa = model[aName, "Mean Sq"]
            else warning(paste("missing factor", aName, "in model"))
            if ("Residuals" %in% row.names(model)) 
                MSe = model["Residuals", "Mean Sq"]
            else warning("missing Residuals in model")
            Cb = Ca = Cab = Cerror = 0
            Cb = (MSb - MSe)/(numA * numMPP)
            Ca = (MSa - MSe)/(numB * numMPP)
            Cab = 0
            Cerror = (MSe)
            gdo@RedANOVA = redFit
        }
        gdo@method = "crossed"
        Ca = max(0, Ca)
        Cb = max(0, Cb)
        Cab = max(0, Cab)
        totalRR = Ca + Cab + Cerror
        repeatability = Cerror
        reproducibility = Ca + Cab
        bTob = max(0, Cb)
        totalVar = Cb + Ca + Cab + Cerror
        estimates = list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
        varcomp = list(totalRR = totalRR, repeatability = repeatability, reproducibility = reproducibility, a = Ca, a_b = Cab, bTob = bTob, totalVar = totalVar)
        gdo@Estimates = estimates
        gdo@Varcomp = varcomp
    }
    cat("\n")
    cat(paste("AnOVa Table - ", gdo@method, "Design\n"))
    print(summary(gdo@ANOVA))
    cat("\n")
    cat("----------\n")
    if (!identical(gdo@RedANOVA, gdo@ANOVA) && gdo@method == "crossed") {
        cat(paste("AnOVa Table Without Interaction - ", gdo@method, "Design\n"))
        print(summary(gdo@RedANOVA))
        cat("\n")
        cat("----------\n")
    }
    Source = names(gdo@Varcomp)
    Source[Source == "repeatability"] = " repeatability"
    Source[Source == "reproducibility"] = " reproducibility"
    Source[Source == "a_b"] = paste("  ", abName)
    Source[Source == "a"] = paste("  ", aName)
    Source[Source == "bTob"] = bTobName
    VarComp = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))]), 3)
    Contribution = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))])/as.numeric(gdo@Varcomp[length(gdo@Varcomp)]), 3)
    VarComp = t(data.frame(gdo@Varcomp))
    VarCompContrib = VarComp/gdo@Varcomp$totalVar
    Stdev = sqrt(VarComp)
    StudyVar = Stdev * gdo@Sigma
    StudyVarContrib = StudyVar/StudyVar["totalVar", ]
    SNR = 1
    ptRatio = NULL
    temp = NULL
    if ((length(gdo@GageTolerance) > 0) && (gdo@GageTolerance > 0)) {
        ptRatio = StudyVar/gdo@GageTolerance
        temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
        names(temp)[6] = c("P/T Ratio")
        row.names(temp) = c(Source)
    }
    else {
        temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
        row.names(temp) = c(Source)
    }
    cat("\n")
    cat("Gage R&R\n")
    tempout = temp
    print(format(tempout, digits = dig))
    cat("\n")
    cat("---\n")
    cat(" * Contrib equals Contribution in %\n")
    SNRTemp = sqrt(2) * (temp[bTobName, "Stdev"]/temp["totalRR", "Stdev"])
    if (SNRTemp > 1) 
        SNR = SNRTemp
    cat(paste(" **Number of Distinct Categories (truncated signal-to-noise-ratio) =", floor(SNR), "\n"))
    cat("\n")
    invisible(gdo)
}
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", signature(x = "gageRR"), function(x, y, main, xlab, ylab, col, lwd, fun = mean, ...) {
    horiz = FALSE
    parList = list(...)
    gdo = x
    yName = names(gdo)[5]
    aName = names(gdo)[3]
    bName = names(gdo)[4]
    abName = paste(aName, ":", bName, sep = "")
    if (missing(col)) 
        col = 2:(length(unique(gdo[, 3])) + 1)
    if (missing(lwd)) 
        lwd = 1
    par(mfrow = c(3, 2))
    temp = NULL
    Source = names(gdo@Varcomp)
    VarComp = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))]), 3)
    Contribution = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))])/as.numeric(gdo@Varcomp[length(gdo@Varcomp)]), 3)
    VarComp = t(data.frame(gdo@Varcomp))
    VarCompContrib = VarComp/gdo@Varcomp$totalVar
    Stdev = sqrt(VarComp)
    StudyVar = Stdev * gdo@Sigma
    StudyVarContrib = StudyVar/StudyVar["totalVar", ]
    if ((length(gdo@GageTolerance) > 0) && (gdo@GageTolerance > 0)) {
        ptRatio = StudyVar/gdo@GageTolerance
        temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib, ptRatio)
        contribFrame = data.frame(VarCompContrib, StudyVarContrib, ptRatio)
        names(temp)[6] = c("P/T Ratio")
        row.names(temp) = c(Source)
        SNR = sqrt(2 * (temp["bTob", "VarComp"]/temp["totalRR", "VarComp"]))
    }
    else {
        temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar, StudyVarContrib)
        contribFrame = data.frame(VarCompContrib, StudyVarContrib)
    }
    bTob = paste(bName, "To", bName, sep = "")
    Source[Source == "bTob"] = bTob
    row.names(contribFrame) = Source
    if (gdo@method == "crossed") 
        contribFrame = contribFrame[-match(c("totalVar", "a", "a_b"), row.names(temp)), ]
    else contribFrame = contribFrame[-match(c("totalVar"), row.names(temp)), ]
    numBars = ncol(contribFrame)
    ymax = max(max(contribFrame))
    main1 = NA
    if (missing(main) || is.na(main[1])) 
        main1 = "Components of Variation"
    else main1 = main[1]
    xlab1 = NA
    if (missing(xlab) || is.na(xlab[1])) 
        xlab1 = "component"
    else xlab1 = xlab[1]
    ylab1 = NA
    if (missing(ylab) || is.na(ylab[1])) 
        ylab1 = ""
    else ylab1 = ylab[1]
    argList = list(...)
    redList = argList[names(argList) != "cex"]
    mybp = do.call(barplot, c(list(t(contribFrame), xlab = xlab1, ylab = ylab1, main = main1, names.arg = rep("", 4), axes = F, beside = T, ylim = c(0, 1.3 * 
        ymax), col = col[1:numBars]), redList))
    axis(1, at = colMeans(mybp), labels = names(as.data.frame(t(contribFrame))), ...)
    axis(2, ...)
    box()
    legend("topright", names(contribFrame), col = col[1:numBars], pch = c(15, 15), horiz = horiz, inset = 0.02)
    if (gdo@method == "crossed") {
        main2 = NA
        if (missing(main) || is.na(main[2])) 
            main2 = paste(yName, "by", bName)
        else main2 = main[2]
        xlab2 = NA
        if (missing(xlab) || is.na(xlab[2])) 
            xlab2 = bName
        else xlab2 = xlab[2]
        ylab2 = NA
        if (missing(ylab) || is.na(ylab[2])) 
            ylab2 = yName
        else ylab2 = ylab[2]
        boxplot(split(gdo[, yName], gdo[, bName]), xlab = xlab2, ylab = ylab2, main = main2, ...)
        mByPa = split(gdo[, 5], as.numeric(gdo[, 4]))
        lines(sort(as.numeric(gdo[, 4])), lapply(mByPa, median)[sort(as.numeric(gdo[, 4]))], lwd = lwd)
        points(sort(as.numeric(gdo[, 4])), lapply(mByPa, median)[sort(as.numeric(gdo[, 4]))], lwd = lwd, pch = 13, cex = 2)
        main3 = NA
        if (missing(main) || is.na(main[3])) 
            main3 = paste(yName, "by", aName)
        else main3 = main[3]
        xlab3 = NA
        if (missing(xlab) || is.na(xlab[3])) 
            xlab3 = aName
        else xlab3 = xlab[3]
        ylab3 = NA
        if (missing(ylab) || is.na(ylab[3])) 
            ylab3 = yName
        else ylab3 = ylab[3]
        colVec = .mapping(gdo[, 3], sort(unique(gdo[, 3])), col[1:length(unique(gdo[, 3]))])
        boxplot(split(gdo[, yName], gdo[, aName]), col = colVec, xlab = xlab3, ylab = ylab3, main = main3, ...)
        mByOp = split(gdo[, 5], as.numeric(gdo[, 3]))
        lines(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd)
        points(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, median)[sort(names(mByOp))], lwd = lwd, pch = 13, cex = 2)
        agg = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = mean)
        tab = table(agg[, 2])
        sgSize = tab[1]
        aggSd = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = sd)
        tab = table(aggSd[, 2])
        sm = mean(aggSd[, 3])
        aggMean = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
        xm = mean(agg[, 3])
        UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
        LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
        values = c(UCL, xm, LCL)
        old.par = par()$mar
        par(mar = c(5.1, 4.1, 4.1, 10.1))
        plot(agg[, 3], type = "n", axes = FALSE, xlab = aName, ylab = expression(bar(x)), main = expression(paste(bar(x), " Chart")))
        box()
        abline(h = xm, col = 3)
        abline(h = UCL, col = 2)
        abline(h = LCL, col = 2)
        axis(2)
        axis(4, at = c(xm, UCL, LCL), labels = c("", "", ""))
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, xm, substitute(bar(x) == xm, list(xm = round(xm, 2))), adj = 0, srt = 0, xpd = TRUE)
        par(mar = old.par)
        j = 1
        for (i in 1:length(tab)) {
            lines(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
            points(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
            if (i < length(tab)) 
                abline(v = j + tab[i] - 1 + 0.5, lty = 2)
            axis(1, at = j, labels = names(tab[i]))
            j = j + tab[i]
        }
        main4 = NA
        if (missing(main) || is.na(main[4])) 
            main4 = paste("Interaction", abName)
        else main4 = main[4]
        xlab4 = NA
        if (missing(xlab) || is.na(xlab[4])) 
            xlab4 = names(gdo)[4]
        else xlab4 = xlab[4]
        ylab4 = NA
        if (missing(ylab) || is.na(ylab[4])) 
            ylab4 = paste(as.character(body(match.fun(fun)))[2], "of", names(gdo)[5])
        else ylab4 = ylab[4]
        old.par = par()$mar
        par(mar = c(5.1, 4.1, 4.1, 10.1))
        .aip(gdo[, 4], gdo[, 3], response = gdo[, 5], xlab = xlab4, ylab = ylab4, main = main4, col = col, type = "b", title = names(gdo)[3], ...)
        par(mar = old.par)
        D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
        D4 = c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
        helpRange = function(x) {
            return(diff(range(x)))
        }
        aggForLimits = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = helpRange)
        Rm = mean(aggForLimits[, 3])
        UCL = D4[sgSize] * Rm
        LCL = D3[sgSize] * Rm
        agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = helpRange)
        tab = table(agg[, 2])
        sgSize = tab[1]
        old.par = par()$mar
        par(mar = c(5.1, 4.1, 4.1, 10.1))
        plot(agg[, 3], ylim = c(0, max(max(agg[, 3]), UCL)), type = "n", xlab = aName, ylab = "R", axes = FALSE, main = "R Chart")
        axis(2)
        axis(4, at = c(Rm, UCL, LCL), labels = c("", "", ""))
        box()
        abline(h = Rm, col = 3)
        abline(h = UCL, col = 2)
        abline(h = LCL, col = 2)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, Rm, substitute(bar(R) == Rm, list(Rm = round(Rm, 2))), adj = 0, srt = 0, xpd = TRUE)
        j = 1
        for (i in 1:length(tab)) {
            lines(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
            points(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
            if (i < length(tab)) 
                abline(v = j + tab[i] - 1 + 0.5, lty = 2)
            axis(1, at = j, labels = names(tab[i]))
            j = j + tab[i]
        }
        par(mar = old.par)
    }
    if(gdo@method == "nested")
    {
        main2 = NA
        if (missing(main) || is.na(main[2])) 
            main2 = paste(yName, "By", bName, "Within", aName)
        else main2 = main[2]
        xlab2 = NA
        if (missing(xlab) || is.na(xlab[2])) 
            xlab2 = NA
        else xlab2 = xlab[2]
        ylab2 = NA
        if (missing(ylab) || is.na(ylab[2])) 
            ylab2 = yName
        else ylab2 = ylab[2]
        agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = function(x) {
            return(x)
        })
        plot(1:nrow(agg), main = main2, xlab = xlab2, ylab = ylab2, ylim = range(agg[, 3]), axes = FALSE)
        axis(2)
        box()
        label2 = ""
        for (i in 1:nrow(agg)) {
            points(rep(i, length(agg[i, 3])), agg[i, 3])
            axis(1, at = i, labels = agg[i, 1])
            if (agg[i, 2] != label2) {
                axis(1, at = i, labels = agg[i, 2], line = 1, tick = FALSE)
                label2 = agg[i, 2]
            }
        }
        aggm = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
        lines(aggm[, 3])
        points(aggm[, 3], pch = 13, cex = 2)
        main3 = NA
        if (missing(main) || is.na(main[3])) 
            main3 = paste(yName, "by", aName)
        else main3 = main[3]
        xlab3 = NA
        if (missing(xlab) || is.na(xlab[3])) 
            xlab3 = aName
        else xlab3 = xlab[3]
        ylab3 = NA
        if (missing(ylab) || is.na(ylab[3])) 
            ylab3 = yName
        else ylab3 = ylab[3]
        colVec = .mapping(gdo[, 3], sort(unique(gdo[, 3])), col[1:length(unique(gdo[, 3]))])
        boxplot(split(gdo[, yName], gdo[, aName]), col = colVec, xlab = xlab3, ylab = ylab3, main = main3, ...)
        mByOp = split(gdo[, 5], as.numeric(gdo[, 3]))
        lines(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd)
        points(sort(as.numeric(factor(names(mByOp)))), lapply(mByOp, mean)[sort(names(mByOp))], lwd = lwd, pch = 13, cex = 2)
        
                agg = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = mean)
        tab = table(agg[, 2])
        sgSize = tab[1]
        aggSd = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = sd)
        tab = table(aggSd[, 2])
        sm = mean(aggSd[, 3])
        aggMean = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = mean)
        xm = mean(agg[, 3])
        UCL = xm + ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
        LCL = xm - ((3 * sm)/(.c4(sgSize) * sqrt(sgSize)))
        values = c(UCL, xm, LCL)
        old.par = par()$mar
        par(mar = c(5.1, 4.1, 4.1, 10.1))
        plot(agg[, 3], type = "n", axes = FALSE, xlab = aName, ylab = expression(bar(x)), main = expression(paste(bar(x), " Chart")))
        box()
        abline(h = xm, col = 3)
        abline(h = UCL, col = 2)
        abline(h = LCL, col = 2)
        axis(2)
        axis(4, at = c(xm, UCL, LCL), labels = c("", "", ""))
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, xm, substitute(bar(x) == xm, list(xm = round(xm, 2))), adj = 0, srt = 0, xpd = TRUE)
        par(mar = old.par)
        j = 1
        for (i in 1:length(tab)) {
            lines(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
            points(j:(j + tab[i] - 1), aggMean[j:(j + tab[i] - 1), 3])
            if (i < length(tab)) 
                abline(v = j + tab[i] - 1 + 0.5, lty = 2)
            axis(1, at = j, labels = names(tab[i]))
            j = j + tab[i]
        }

                par(mar = old.par)
        D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.284, 0.308, 0.329, 0.348)
        D4 = c(0, 3.267, 2.574, 2.282, 2.115, 2.004, 1.924, 1.864, 1.816, 1.777, 1.744, 1.716, 1.692, 1.671, 1.652)
        helpRange = function(x) {
            return(diff(range(x)))
        }
        aggForLimits = aggregate(gdo[, yName], list(gdo[, aName], gdo[, bName]), FUN = helpRange)
        Rm = mean(aggForLimits[, 3])
        UCL = D4[sgSize] * Rm
        LCL = D3[sgSize] * Rm
        agg = aggregate(gdo[, yName], list(gdo[, bName], gdo[, aName]), FUN = helpRange)
        tab = table(agg[, 2])
        sgSize = tab[1]
        old.par = par()$mar
        par(mar = c(5.1, 4.1, 4.1, 10.1))
        plot(agg[, 3], ylim = c(0, max(max(agg[, 3]), UCL)), type = "n", xlab = aName, ylab = "R", axes = FALSE, main = "R Chart")
        axis(2)
        axis(4, at = c(Rm, UCL, LCL), labels = c("", "", ""))
        box()
        abline(h = Rm, col = 3)
        abline(h = UCL, col = 2)
        abline(h = LCL, col = 2)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, LCL, paste("LCL =", round(LCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, UCL, paste("UCL =", round(UCL, 2)), adj = 0, srt = 0, xpd = TRUE)
        text(length(agg[, 3]) + length(agg[, 3]) * 0.075, Rm, substitute(bar(R) == Rm, list(Rm = round(Rm, 2))), adj = 0, srt = 0, xpd = TRUE)
        j = 1
        for (i in 1:length(tab)) {
            lines(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
            points(j:(j + tab[i] - 1), agg[j:(j + tab[i] - 1), 3])
            if (i < length(tab)) 
                abline(v = j + tab[i] - 1 + 0.5, lty = 2)
            axis(1, at = j, labels = names(tab[i]))
            j = j + tab[i]
        
         }
        
        
        
#        plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
#        plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
        plot(1, 1, type = "n", axes = FALSE, xlab = NA, ylab = NA, main = NA)
    }
}) 
