setClass(Class = "mixDesign", representation = representation(name = "character", factors = "list", total = "numeric", lower = "numeric", design = "data.frame", 
    designType = "character", pseudo = "data.frame", response = "data.frame", Type = "data.frame", block = "data.frame", runOrder = "data.frame", standardOrder = "data.frame", 
    desireVal = "list", desirability = "list", fits = "data.frame"))
setMethod("factors", "mixDesign", function(x) x@factors)
setReplaceMethod("factors", "mixDesign", function(x, value) {
    if (length(value) != ncol(x@pseudo)) 
        stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
    x@factors <- value
    x
})
setMethod("names", "mixDesign", function(x) {
    return(sapply(x@factors, names))
})
setReplaceMethod("names", "mixDesign", function(x, value) {
    for (i in 1:length(x@factors)) names(x@factors[[i]]) = as.character(value[i])
    x
})
setMethod("as.data.frame", "mixDesign", function(x, row.names = NULL, optional = FALSE) {
    frameOut = cbind(x@standardOrder, x@runOrder, x@Type, x@pseudo, x@response)
    return(frameOut)
})
as.data.frame.mixDesign = function(x, row.names = NULL, optional = FALSE, ...) {
    frameOut = cbind(x@standardOrder, x@runOrder, x@Type, x@pseudo, x@response)
    return(frameOut)
}
setMethod("show", signature(object = "mixDesign"), function(object) {
    print(format(as.data.frame(object), digits = 4))
})
setMethod("response", "mixDesign", function(object) {
    return(object@response)
})
setReplaceMethod("response", "mixDesign", function(object, value) {
    print(deparse(substitute(value)))
    if (!is.numeric(value) & !is.data.frame(value)) 
        stop("vector or data.frame must be given")
    if (is.numeric(value)) {
        if (length(value) != nrow(object@pseudo)) 
            stop("differing lengths")
        temp = data.frame(value)
        names(temp) = deparse(substitute(value))[1]
        value = temp
    }
    if (is.data.frame(value)) {
        if (nrow(value) != nrow(object@pseudo)) 
            stop("differing number of rows")
    }
    object@response = value
    return(object)
})
.npp = function(mdo) {
    pseudo = mdo@pseudo
    Type = mdo@Type
    temp = as.character(as.data.frame(mdo)$Type)
    tab = table(temp)
    nums = data.frame(matrix(tab, nrow = 4, ncol = length(tab), byrow = TRUE))
    nums[1:2, ] = NA
    nums[4, ] = c(nrow(Type), rep("", length(tab) - 1))
    row.names(nums) = c("Unique", "Replicates", "Sub Total", "Total")
    names(nums) = names(tab)
    for (i in names(tab)) {
        sSet = pseudo[Type == i, ]
        usSet = unique(sSet)
        nums["Unique", i] = nrow(usSet)
        if (nrow(usSet) == 1) {
            nums["Replicates", i] = nrow(sSet)
        }
        else {
            for (j in 1:nrow(usSet)) {
                uCount1 = sum(apply(apply(sSet, 1, "==", usSet[j, ]), 2, "all") * 1)
                if (j == 1) {
                  uCount2 = uCount1
                  nums["Replicates", i] = uCount1
                }
                if (uCount2 != uCount1) 
                  nums["Replicates", i] = -1
            }
        }
    }
    cat("Information about the Design Points:")
    cat("\n")
    cat("\n")
    print(nums)
}
setMethod("show", signature(object = "mixDesign"), function(object) {
    print(as.data.frame(object))
})
setMethod(".nfp", "mixDesign", function(object) {
    x = factors(object)
    if (is.list(x) && length(x[[1]]) > 0) {
        numAttr = length(attributes(x[[1]]))
        .numFac = length(x)
        frameOut = data.frame(matrix(ncol = .numFac, nrow = (numAttr - 1)))
        for (i in 1:(numAttr - 1)) {
            charVec = character(0)
            for (j in 1:.numFac) {
                charVec = c(charVec, names(attributes(x[[1]])[i]), "\t\t")
                frameOut[i, j] = attributes(x[[j]])[[i]]
            }
        }
        names(frameOut) = names(x)
        rownames(frameOut) = names(attributes(x[[1]]))[1:(numAttr - 1)]
    }
    else {
        stop("no list given or length of list < 1")
    }
    print(frameOut)
})
setMethod("summary", signature(object = "mixDesign"), function(object) {
    cat(paste("Simplex", toupper(object@designType), "Design"))
    cat("\n")
    cat("Information about the factors:\n\n")
    .nfp(object)
    cat("\n-----------\n")
    cat("\n")
    .npp(object)
    cat("\n-----------\n")
    cat("\n")
    cat("Information about the constraints:\n\n")
    lower = object@lower
    temp = character(0)
    for (i in seq(along = lower)) temp = c(temp, paste(LETTERS[i], ">=", lower[i]))
    cat(temp)
    cat("\n")
    cat("\n-----------\n")
    cat("\n")
    times = nrow(object@pseudo)
    pseudo = format(object@pseudo, digits = 2)
    design = format(object@design, digits = 2)
    amount = design
    if (object@total[2] != 1) 
        amount = format(object@design * object@total[2], digits = 2)
    temp = c("                             ", "PseudoComponent", "_|_", "Proportion", "_|_", "Amount")
    cat(temp)
    cat("\n")
    cat("\n")
    temp = cbind(pseudo, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), design)
    temp = cbind(temp, `_` = rep(" ", times = times), `|` = rep("|", times = times), `_` = rep(" ", times = times), amount)
    temp = cbind(object@standardOrder, object@runOrder, object@Type, `|` = rep("|", times = times), temp, `|` = rep("|", times = times), object@response)
    show(temp)
    cat("\n-----------\n")
    cat("\n")
    cat(paste("Mixture Total:", object@total[1], "equals", object@total[2]))
    cat("\n")
    cat("\n")
    invisible(as.data.frame(object))
})
setMethod("units", "mixDesign", function(x) {
    return(sapply(factors(x), .unit))
})
setMethod("units<-", "mixDesign", function(x, value) {
    for (i in 1:length(x@factors)) if (length(value) > 1) 
        .unit(x@factors[[i]]) = as.character(value[i])
    else .unit(x@factors[[i]]) = as.character(value[1])
    x
})
setMethod("highs", "mixDesign", function(object) {
    listOut = vector(mode = "list")
    for (i in names(factors(object))) {
        listOut[i] = .high(object@factors[[i]])
    }
    return(listOut)
})
setReplaceMethod("highs", "mixDesign", function(object, value) {
    for (i in seq(along = object@factors)) if (length(value) > 1) 
        .high(object@factors[[i]]) = value[i]
    else .high(object@factors[[i]]) = value[1]
    return(object)
})
setMethod("lows", "mixDesign", function(object) {
    listOut = vector(mode = "list")
    for (i in names(factors(object))) {
        listOut[i] = .low(object@factors[[i]])
    }
    return(listOut)
})
setReplaceMethod("lows", "mixDesign", function(object, value) {
    for (i in seq(along = object@factors)) {
        if (length(value) > 1) 
            .low(object@factors[[i]]) = value[i]
        else .low(object@factors[[i]]) = value[1]
    }
    return(object)
})
.permut = function(x) {
    if (any(is.na(x))) 
        stop(paste(deparse(substitute(x)), "contains NA"))
    x = sort(x, decreasing = FALSE)
    n = length(x)
    num = 1:n
    frameOut = matrix(NA, nrow = 1, ncol = n)
    frameOut[1, ] = x
    while (TRUE) {
        highest = NA
        for (j in 1:(n - 1)) {
            if (x[j] < x[j + 1]) 
                highest = j
        }
        if (is.na(highest)) 
            return(frameOut)
        else {
            l = max((num)[x[highest] < x])
            temp = x[l]
            x[l] = x[highest]
            x[highest] = temp
            x[(highest + 1):n] = rev(x[(highest + 1):n])
        }
        frameOut = rbind(frameOut, x, deparse.level = 2)
    }
}
.simplexCentroid = function(p) {
    if (p <= 1 | !is.numeric(p)) 
        stop("invalid value for p")
    frameOut = NA
    for (i in 1:p) {
        initial = rep(0, times = p)
        initial[1:i] = 1/i
        mat = .permut(initial)
        if (i == 1) 
            frameOut = mat
        else frameOut = rbind(frameOut, mat, deparse.level = 2)
    }
    frameOut = data.frame(frameOut, row.names = 1:(2^p - 1))
    names(frameOut) = LETTERS[1:p]
    return(frameOut)
}
.simplexCentroid(4)
mixDesign = function(p, n = 3, type = "lattice", center = TRUE, axial = FALSE, delta, replicates = 1, lower, total = 1, randomize, seed) {
    DB = FALSE
    frameOut = NA
    out = new("mixDesign")
    if (missing(p)) 
        stop("the number of factors p must be given")
    if (p <= 1 | !is.numeric(p)) 
        stop("invalid value for p")
    if (!(type %in% c("lattice", "centroid"))) 
        stop("type needs to be \"lattice\" or \"centroid\"")
    out@designType = type
    if (missing(delta)) 
        delta = (p - 1)/(2 * p)
    if (missing(randomize)) 
        randomize = TRUE
    if (!missing(seed)) 
        set.seed(seed)
    if (!is.numeric(total)) 
        stop("total needs to be a numeric vector with <= 2 arguments")
    else {
        if (total[1] > 1 || total[1] <= 0) 
            stop("totall[1] needs to be within (0,1]")
        if (is.na(total[2])) 
            total[2] = 1
        if (total[2] <= 0) 
            stop("total[2] needs to be > 0")
    }
    out@total = total
    if (!is.numeric(replicates)) 
        stop("replicates need to be numeric")
    if (delta > (p - 1)/p) {
        delta = (p - 1)/(2 * p)
        warning(paste("delta was reset to:", delta))
    }
    if (missing(lower)) 
        lower = 0
    if (length(lower) == 1) 
        lower = rep(lower, p)
    if (length(lower) > 1) {
        initial = rep(0, p)
        if (length(lower) < p) 
            lower[, (length(lower) + 1):p] = 0
    }
    out@lower = lower
    repTemp = list(center = 1, axial = 1)
    for (i in 1:(p - 1)) repTemp[[as.character(i)]] = 1
    if (length(replicates) > 1) {
        for (i in 1:length(replicates)) repTemp[[i]] = replicates[[i]]
        replicates = repTemp
    }
    if (length(replicates) == 1) {
        for (i in 1:length(repTemp)) repTemp[[i]] = replicates
        replicates = repTemp
    }
    if (DB) 
        print(replicates)
    N = factorial(p + n - 1)/(factorial(n) * factorial(p - 1))
    if (identical(type, "lattice")) {
        j = 1
        x = numeric(p)
        j = 1
        x[1] = n
        for (i in 1:N) {
            x[j + 1] = n - (sum(x[1:j]))
            if (j < (p - 1)) 
                x[(j + 2):p] = 0
            if (i == 1) 
                frameOut = data.frame(matrix(x, ncol = p, nrow = 1))
            else frameOut = rbind(frameOut, x)
            if (DB) 
                print(x)
            logVec = rep(FALSE, p)
            logVec[1:(p - 1)] = x[1:p - 1] > 0
            if (any(logVec)) {
                j = max((1:p)[logVec])
                x[j] = x[j] - 1
            }
        }
        frameOut = (frameOut/(n))
        names(frameOut) = LETTERS[1:p]
    }
    if (identical(type, "centroid")) {
        frameOut = .simplexCentroid(p)
    }
    frameOutCopy = frameOut[1, ]
    Type = data.frame(Type = "1")
    temp = apply(ceiling(frameOut), 2, ">", 0) * 1
    temp = apply(temp, 1, sum)
    for (i in 1:nrow(frameOut)) {
        typ = as.character(temp[i])
        times = replicates[[typ]]
        if (is.null(times)) 
            times = 0
        else times = times - 1
        repFrame = frameOut[i, ]
        if (all(frameOut[i, ] == frameOut[i, 1])) 
            typFrame = data.frame(Type = "center")
        else typFrame = data.frame(Type = paste(typ, "-blend", sep = ""))
        if (times >= 1) {
            for (j in 1:times) {
                repFrame = rbind(repFrame, frameOut[i, ])
                typFrame = rbind(typFrame, data.frame(Type = paste(typ, "-blend", sep = "")))
            }
        }
        frameOutCopy = rbind(frameOutCopy, repFrame)
        Type = rbind(Type, typFrame)
        if (i == 1) {
            frameOutCopy = frameOutCopy[c(-1), ]
            Type = data.frame(Type = Type[c(-1), ])
        }
    }
    frameOut = frameOutCopy
    if (DB) {
        print(Type)
        print(frameOutCopy)
    }
    keepIndex = (1:nrow(frameOut))[!apply(Type, 1, "==", "center")]
    Type = data.frame(Type = Type[keepIndex, ])
    frameOut = frameOut[keepIndex, ]
    if (center) {
        center = data.frame(matrix(1/p, nrow = 1, ncol = p))
        if (DB) 
            print(center)
        times = replicates$center
        if (is.null(times)) 
            times = 0
        else times = times - 1
        if (n == p) 
            times = times - 1
        if ((n%%p) == 0) {
            numCenter = n%/%p
        }
        temp = center
        if (times >= 1) {
            for (i in 1:times) center = rbind(center, temp)
        }
        names(center) = names(frameOut)
        frameOut = rbind(frameOut, center)
        Type = rbind(Type, data.frame(Type = rep("center", times + 1)))
        if (DB) 
            print(frameOut)
    }
    if (axial) {
        temp = rep(NA, p)
        axial = data.frame(matrix(NA, ncol = p, nrow = p))
        for (i in 1:p) {
            temp[i] = delta + 1/p
            temp[c(-i)] = (1 - temp[i])/(p - 1)
            axial[i, ] = temp
        }
        times = replicates$axial
        if (is.null(times)) 
            times = 0
        else times = times - 1
        temp = axial
        if (times >= 1) {
            for (i in 1:times) axial = rbind(axial, temp)
        }
        names(axial) = names(frameOut)
        frameOut = rbind(frameOut, axial)
        Type = rbind(Type, data.frame(Type = rep("axial", (times + 1) * p)))
        if (DB) {
            print(frameOut)
            print(Type)
        }
    }
    StandOrder = 1:nrow(frameOut)
    RunOrder = StandOrder
    if (randomize) {
        RunOrder = sample(1:nrow(frameOut), nrow(frameOut), replace = FALSE, prob = NULL)
    }
    frameOut = frameOut[order(RunOrder), ]
    row.names(frameOut) = frameOut$RunOrder
    out@pseudo = frameOut
    out@runOrder = data.frame(RunOrder = RunOrder)
    out@standardOrder = data.frame(StandOrder = StandOrder)
    out@Type = data.frame(Type = Type[order(RunOrder), ])
    out@response = data.frame(y = rep(NA, nrow(out@pseudo)))
    design = frameOut
    design[, ] = NA
    for (i in 1:ncol(frameOut)) {
        design[, i] = frameOut[, i] * (total[1] - sum(lower)) + lower[i]
    }
    out@design = design
    listFac = vector("list", p)
    for (i in seq(along = listFac)) listFac[i] = new("doeFactor")
    names(listFac) = LETTERS[1:p]
    factors(out) = listFac
    if (out@total[2] != 1) {
        lows(out) = lower * out@total[2]
        highs(out) = 1 * out@total[2]
        units(out) = "NA"
    }
    else if (any(out@lower != 0)) {
        lows(out) = out@lower
        highs(out) = 1 * out@total[1]
        units(out) = "%"
    }
    else {
        lows(out) = 0
        highs(out) = 1 * out@total[1]
        units(out) = "%"
    }
    return(out)
} 
