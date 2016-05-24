.NAMES = LETTERS[c(1:8, 10:26)]
setClass("doeFactor", representation = representation(low = "ANY", high = "ANY", name = "character",
    unit = "character", type = "character"), prototype = prototype(low = -1, high = 1, name = "",
    unit = "", type = "numeric"))
setGeneric(".low", function(object) standardGeneric(".low"))
setGeneric(".low<-", function(x, value) standardGeneric(".low<-"))
setMethod(".low", "doeFactor", function(object) unlist(object@low))
setReplaceMethod(".low", "doeFactor", function(x, value) {
    boolOld = is.numeric(.low(x))
    x@low <- value
    boolNew = is.numeric(.low(x))
    if (boolNew)
        x@type = "numeric"
    else x@type = "factor"
    if (boolOld != boolNew)
        print("Note: The types of the factors were changed!")
    return(x)
})
setGeneric(".high", function(object) standardGeneric(".high"))
setGeneric(".high<-", function(x, value) standardGeneric(".high<-"))
setMethod(".high", "doeFactor", function(object) unlist(object@high))
setReplaceMethod(".high", "doeFactor", function(x, value) {
    boolOld = is.numeric(.high(x))
    x@high <- value
    boolNew = is.numeric(.high(x))
    if (boolNew)
        x@type = "numeric"
    else x@type = "factor"
    if (boolOld != boolNew)
        print("Note: The types of the factors were changed!")
    return(x)
})
code2real = function(low, high, codedValue) {
    return((diff(c(low, high))/2) * codedValue + mean(c(low, high)))
}
code2real(10, 20, -1)
setGeneric(".type", function(object) standardGeneric(".type"))
setGeneric(".type<-", function(x, value) standardGeneric(".type<-"))
setMethod(".type", "doeFactor", function(object) object@type)
setReplaceMethod(".type", "doeFactor", function(x, value) {
    x@type <- value
    x
})
setGeneric(".unit", function(object) standardGeneric(".unit"))
setGeneric(".unit<-", function(x, value) standardGeneric(".unit<-"))
setMethod(".unit", "doeFactor", function(object) object@unit)
setReplaceMethod(".unit", "doeFactor", function(x, value) {
    x@unit <- value
    x
})
setMethod("names", "doeFactor", function(x) {
    x@name
})
setReplaceMethod("names", "doeFactor", function(x, value) {
    x@name <- value
    x
})
setMethod("show", signature(object = "doeFactor"), function(object) {
    cat("Name: ", names(object), "\n")
    cat("low Setting: ", .low(object), "\n")
    cat("high setting: ", .high(object), "\n")
    cat("Unit: ", .unit(object), "\n")
    cat("type: ", .type(object), "\n")
    cat("\n")
})
setClass(Class = "facDesign", representation = representation(name = "character", factors = "list",
    cube = "data.frame", star = "data.frame", centerCube = "data.frame", centerStar = "data.frame",
    generator = "ANY", response = "data.frame", block = "data.frame", blockGen = "data.frame", runOrder = "data.frame",
    standardOrder = "data.frame", desireVal = "list", desirability = "list", fits = "list"))
setGeneric(".nfp", function(object) standardGeneric(".nfp"))
setMethod(".nfp", "facDesign", function(object) {
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
setGeneric("fits", function(x) standardGeneric("fits"))
setGeneric("fits<-", function(x, value) standardGeneric("fits<-"))
setMethod("fits", "facDesign", function(x) {
    return(x@fits)
})
setMethod("fits<-", "facDesign", function(x, value) {
    if (!identical(class(value), "lm"))
        stop(paste(deparse(substitute(value)), "needs to an object of class lm"))
    if (!any(names(value$model)[1] == names(response(x))))
        stop(paste("fitted response", names(value$model)[1], "could not be found in", deparse(substitute(x))))
    listPos = length(x@fits) + 1
    yName = names(value$model)[1]
    isIn = (yName == names(x@fits))
    if (any(isIn))
        listPos = (1:length(names(x@fits)))[isIn]
    x@fits[[listPos]] = value
    names(x@fits)[listPos] = yName
    x
})
setGeneric("desires", function(x) standardGeneric("desires"))
setGeneric("desires<-", function(x, value) standardGeneric("desires<-"))
setMethod("desires", "facDesign", function(x) {
    return(x@desirability)
})
setMethod("desires<-", "facDesign", function(x, value) {
    if (!any(value@response == names(response(x))))
        stop(paste(value@response, "is not a response!"))
    listPos = length(x@desirability) + 1
    yName = value@response
    isIn = (yName == names(x@desirability))
    if (any(isIn))
        listPos = (1:length(names(x@desirability)))[isIn]
    x@desirability[[listPos]] = value
    names(x@desirability)[listPos] = yName
    x
})
setMethod("nrow", "facDesign", function(x) nrow(as.data.frame(x)))
setMethod("ncol", "facDesign", function(x) ncol(as.data.frame(x)))
.numFac = function(fdo) {
    return(length(names(fdo)))
}
setGeneric(".clear", function(x) standardGeneric(".clear"))
setMethod(".clear", "facDesign", function(x) {
    x@standardOrder = data.frame()
    x@runOrder = data.frame()
    x@cube = data.frame()
    x@centerStar = data.frame()
    x@centerCube = data.frame()
    x@star = data.frame()
    x@block = data.frame()
    x@blockGen = data.frame()
    x@response = data.frame()
    return(x)
})
setMethod("[", signature(x = "facDesign", i = "ANY", j = "ANY"), function(x, i, j) {
    return(as.data.frame(x)[i, j])
})
.helpSort = function(fdo, runOrd = TRUE) {
    oldIndex = standOrd(fdo)
    newIndex = order(runOrd(fdo))
    out = fdo
    cube(out) = oldIndex[newIndex[1:nrow(cube)]]
}
setReplaceMethod("[", signature(x = "facDesign", i = "ANY", j = "ANY"), function(x, i,
    j, value) {
    print(paste("setting values via [", i, ",", j, " ]  is not supported"))
    print("Use response() to set the values of the response!")
    return(x)
    if (class(value) == "data.frame")
        if (dim(as.data.frame(x)) == dim(as.data.frame(value))) {
            print("data.frame in einzelne Slots umkopieren")
        }
    if (!(missing(i) & missing(j))) {
        print("not supported")
    }
    if (missing(j)) {
        print("ncol must be the same")
        return(x)
    }
    if (missing(i)) {
        print("nrow must be the same")
        return(x)
    }
    return(x)
})
setGeneric("types", function(x) standardGeneric("types"))
setGeneric("types<-", function(x, value) standardGeneric("types<-"))
setMethod("types", "facDesign", function(x) {
    return(sapply(factors(x), .type))
})
setReplaceMethod("types", "facDesign", function(x, value) {
    for (i in 1:length(x@factors)) {
        if (!identical(value[i], "numeric") & !identical(value[i], "factor"))
            stop(paste(value[i], "\ttype of factor needs to be 'numeric' or 'factor'"))
        .type(x@factors[[i]]) = as.character(value[i])
    }
    x
})
setMethod("units", "facDesign", function(x) {
    return(sapply(factors(x), .unit))
})
setMethod("units<-", "facDesign", function(x, value) {
    for (i in 1:length(x@factors)) if (length(value) > 1)
        .unit(x@factors[[i]]) = as.character(value[i])
    else .unit(x@factors[[i]]) = as.character(value[1])
    x
})
setGeneric("highs", function(object) standardGeneric("highs"))
setGeneric("highs<-", function(object, value) standardGeneric("highs<-"))
setMethod("highs", "facDesign", function(object) {
    listOut = vector(mode = "list")
    for (i in names(factors(object))) {
        listOut[i] = .high(object@factors[[i]])
    }
    return(listOut)
})
setReplaceMethod("highs", "facDesign", function(object, value) {
    for (i in seq(along = object@factors)) {
        .high(object@factors[[i]]) = value[i]
    }
    return(object)
})
setGeneric("lows", function(object) standardGeneric("lows"))
setGeneric("lows<-", function(object, value) standardGeneric("lows<-"))
setMethod("lows", "facDesign", function(object) {
    listOut = vector(mode = "list")
    for (i in names(factors(object))) {
        listOut[i] = .low(object@factors[[i]])
    }
    return(listOut)
})
setReplaceMethod("lows", "facDesign", function(object, value) {
    for (i in seq(along = object@factors)) {
        .low(object@factors[[i]]) = value[i]
    }
    return(object)
})
setGeneric("cube", function(x) standardGeneric("cube"))
setMethod("cube", "facDesign", function(x) {
    out = x@cube
    out
})
setGeneric("cube<-", function(x, value) standardGeneric("cube<-"))
setReplaceMethod("cube", "facDesign", function(x, value) {
    x@cube <- value
    x
})
setGeneric("star", function(x) standardGeneric("star"))
setGeneric("star<-", function(x, value) standardGeneric("star<-"))
setMethod("star", "facDesign", function(x) x@star)
setReplaceMethod("star", "facDesign", function(x, value) {
    DB = FALSE
    if (!is.data.frame(value))
        stop("data.frame must be provided!")
    if (.numFac(x) != ncol(value))
        stop("number of columns not matching!")
    if (nrow(value) == 0) {
        return("TODO: remove star und Rest anpassen")
    }
    oldResponse = response(x)
    newDf = value
    oldDf = x@star
    numNewRow = nrow(newDf) - nrow(oldDf)
    oldOrd = standOrd(x)
    oldRunOrd = runOrd(x)
    len = nrow(oldOrd)
    lenFirst = nrow(cube(x)) + nrow(centerCube(x))
    standOrd(x) = data.frame(StandOrd = 1:(len + numNewRow))
    newRunOrd = data.frame()
    if (numNewRow > 0) {
        newNums = data.frame(newNums = seq(max(oldRunOrd) + 1, max(oldRunOrd) + numNewRow, by = 1))
        if (DB)
            print(newNums)
        names(newNums) = names(oldRunOrd)
        newRunOrd = data.frame(oldRunOrd[1:lenFirst, ])
        if (DB)
            print(newRunOrd)
        names(newRunOrd) = names(oldRunOrd)
        restFrame = data.frame(oldRunOrd[-c(1:lenFirst), ])
        names(restFrame) = names(oldRunOrd)
        newRunOrd = rbind(newRunOrd, newNums, restFrame)
        if (DB)
            print(newRunOrd)
    }
    else {
        newRunOrd = data.frame(oldRunOrd[1:(lenFirst + nrow(newDf) + nrow(centerStar(x))), ])
        names(newRunOrd) = names(oldRunOrd)
    }
    runOrd(x) = newRunOrd
    naFrame = as.data.frame(matrix(rep(NA, times = ncol(oldResponse) * nrow(newDf)), ncol = ncol(oldResponse)))
    names(naFrame) = names(oldResponse)
    newResponse = data.frame(oldResponse[1:lenFirst, ])
    names(newResponse) = names(oldResponse)
    restFrame = data.frame(oldResponse[-c(1:(lenFirst + nrow(oldDf))), ])
    names(restFrame) = names(oldResponse)
    newResponse = rbind(newResponse, naFrame, restFrame)
    response(x) = newResponse
    if (DB) {
        print(newResponse)
        print("hinter response")
    }
    oldBlockGen = blockGen(x)
    if (ncol(oldBlockGen) > 0) {
        if (DB)
            print("TODO: BlockGen anpassen!")
        newBlockGen = data.frame(oldBlockGen[1:lenFirst, ])
        names(newBlockGen) = names(blockGen(x))
        naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(blockGen(x)) * nrow(newDf)), ncol = ncol(blockGen(x))))
        names(naFrameGen) = names(oldBlockGen)
        restBlockGen = data.frame(oldBlockGen[-c(1:(lenFirst + nrow(oldDf))), ])
        names(restBlockGen) = names(oldBlockGen)
        newBlockGen = rbind(newBlockGen, naFrameGen, restBlockGen)
        if (DB)
            print(newBlockGen)
        blockGen(x) = newBlockGen
    }
    oldBlock = block(x)
    newBlock = data.frame(oldBlock[1:lenFirst, ])
    names(newBlock) = names(oldBlock)
    naFrame = as.data.frame(matrix(rep(max(newBlock) + 1, times = ncol(oldBlock) * nrow(newDf)),
        ncol = ncol(oldBlock)))
    names(naFrame) = names(oldBlock)
    restBlock = data.frame(oldBlock[-c(1:(lenFirst + nrow(oldDf))), ])
    names(restBlock) = names(oldBlock)
    newBlock = rbind(newBlock, naFrame, restBlock)
    block(x) = newBlock
    x@star <- newDf
    return(x)
})
setGeneric("centerCube", function(x) standardGeneric("centerCube"))
setGeneric("centerCube<-", function(x, value) standardGeneric("centerCube<-"))
setMethod("centerCube", "facDesign", function(x) x@centerCube)
setReplaceMethod("centerCube", "facDesign", function(x, value) {
    DB = FALSE
    if (!is.data.frame(value))
        stop("data.frame must be provided!")
    if (.numFac(x) != ncol(value))
        stop("number of columns not matching!")
    if (nrow(value) == 0) {
        return("TODO: remove CenterCube und Rest anpassen")
    }
    newDf = value
    lenCube = nrow(cube(x))
    oldDf = x@centerCube
    oldRunOrd = runOrd(x)
    oldResponse = response(x)
    blockValues = unique(block(x)[1:nrow(cube(x)), ])
    numBlocks = length(blockValues)
    if (numBlocks > 1)
        for (i in 1:(numBlocks - 1)) {
            newDf = rbind(newDf, value)
        }
    if (DB)
        print(newDf)
    numNewRow = nrow(newDf) - nrow(oldDf)
    oldOrd = standOrd(x)
    len = nrow(oldOrd)
    standOrd(x) = data.frame(StandOrd = 1:(len + numNewRow))
    newRunOrd = data.frame()
    if (numNewRow > 0) {
        newNums = data.frame(newNums = seq(max(oldRunOrd) + 1, max(oldRunOrd) + numNewRow, by = 1))
        names(newNums) = names(oldRunOrd)
        if (DB) {
            print("----")
            print(newNums)
        }
        newRunOrd = data.frame(oldRunOrd[1:lenCube, ])
        names(newRunOrd) = names(oldRunOrd)
        restRunOrd = data.frame(oldRunOrd[-c(1:lenCube), ])
        names(restRunOrd) = names(oldRunOrd)
        newRunOrd = rbind(newRunOrd, newNums, restRunOrd)
        if (DB) {
            print("----")
            print(oldRunOrd[-c(1:lenCube), ])
            print("----")
            print(newRunOrd)
        }
        runOrd(x) = newRunOrd
    }
    else {
        newRunOrd = data.frame(oldRunOrd[1:(lenCube + nrow(newDf)), ])
        names(newRunOrd) = names(oldRunOrd)
        restRunOrd = data.frame(oldRunOrd[-c(1:(lenCube + nrow(oldDf))), ])
        names(restRunOrd) = names(oldRunOrd)
        newRunOrd = rbind(newRunOrd, restRunOrd)
        if (DB) {
            print("----")
            print(newRunOrd)
        }
        runOrd(x) = newRunOrd
    }
    naFrame = as.data.frame(matrix(rep(NA, times = ncol(oldResponse) * nrow(newDf)), ncol = ncol(oldResponse)))
    names(naFrame) = names(oldResponse)
    newResponse = data.frame(oldResponse[1:lenCube, ])
    names(newResponse) = names(oldResponse)
    restResponse = data.frame(oldResponse[-c(1:(lenCube + nrow(oldDf))), ])
    names(restResponse) = names(oldResponse)
    newResponse = rbind(newResponse, naFrame, restResponse)
    if (DB) {
        print("newResponse_____")
        print(newResponse)
    }
    response(x) = newResponse
    oldBlockGen = blockGen(x)
    if (ncol(oldBlockGen) > 0) {
        if (DB)
            print("TODO: BlockGen Spalte(n) anpassen")
        newBlockGen = data.frame(oldBlockGen[1:lenCube, ])
        names(newBlockGen) = names(blockGen(x))
        naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(blockGen(x)) * nrow(newDf)), ncol = ncol(blockGen(x))))
        names(naFrameGen) = names(oldBlockGen)
        restFrame = data.frame(oldBlockGen[-c(1:(lenCube + nrow(oldDf))), ])
        names(restFrame) = names(blockGen(x))
        newBlockGen = rbind(newBlockGen, naFrameGen, restFrame)
        blockGen(x) = newBlockGen
        if (DB)
            print(newBlockGen)
    }
    oldBlock = block(x)
    newBlock = data.frame(oldBlock[1:lenCube, ])
    names(newBlock) = names(block(x))
    naFrame = as.data.frame(matrix(rep(blockValues, times = nrow(newDf)/numBlocks), ncol = 1))
    restFrame = as.data.frame(oldBlock[-c(1:(lenCube + nrow(oldDf))), ])
    names(restFrame) = names(block(x))
    if (DB)
        print(naFrame)
    names(naFrame) = names(oldBlock)
    newBlock = rbind(newBlock, naFrame, restFrame)
    block(x) = newBlock
    x@centerCube <- newDf
    return(x)
})
setGeneric("centerStar", function(x) standardGeneric("centerStar"))
setGeneric("centerStar<-", function(x, value) standardGeneric("centerStar<-"))
setMethod("centerStar", "facDesign", function(x) {
    x@centerStar
})
setReplaceMethod("centerStar", "facDesign", function(x, value) {
    DB = FALSE
    if (!is.data.frame(value))
        stop("data.frame must be provided!")
    if (.numFac(x) != ncol(value))
        stop("number of columns not matching!")
    if (nrow(value) == 0) {
        return("TODO: remove CenterCube und Rest anpassen")
    }
    newDf = value
    oldDf = x@centerStar
    numNewRow = nrow(newDf) - nrow(oldDf)
    oldResponse = response(x)
    lenRest = nrow(cube(x)) + nrow(centerCube(x)) + nrow(star(x))
    oldRunOrd = runOrd(x)
    oldOrd = standOrd(x)
    len = nrow(oldOrd)
    standOrd(x) = data.frame(StandOrd = 1:(len + numNewRow))
    newRunOrd = data.frame(oldRunOrd[1:lenRest, ])
    names(newRunOrd) = names(oldRunOrd)
    if (numNewRow > 0) {
        newNums = data.frame(newNums = seq(max(oldRunOrd) + 1, max(oldRunOrd) + numNewRow, by = 1))
        names(newNums) = names(oldRunOrd)
        restFrame = data.frame(oldRunOrd[-c(1:lenRest), ])
        names(restFrame) = names(oldRunOrd)
        newRunOrd = rbind(newRunOrd, newNums, restFrame)
        if (DB)
            print(newRunOrd)
        runOrd(x) = newRunOrd
    }
    else {
        newRunOrd = data.frame(oldRunOrd[1:(lenRest + nrow(newDf)), ])
        names(newRunOrd) = names(oldRunOrd)
        runOrd(x) = newRunOrd
    }
    naFrame = as.data.frame(matrix(rep(NA, times = ncol(oldResponse) * nrow(newDf)), ncol = ncol(oldResponse)))
    names(naFrame) = names(oldResponse)
    newResponse = data.frame(oldResponse[1:lenRest, ])
    names(newResponse) = names(response(x))
    newResponse = rbind(newResponse, naFrame)
    if (DB)
        print(" vor centerStar response")
    response(x) = newResponse
    if (DB)
        print("hinter centerStar response")
    oldBlockGen = blockGen(x)
    if (ncol(oldBlockGen) > 0) {
        print("TODO: BlockGen Spalte(n) anpassen")
        newBlockGen = data.frame(oldBlockGen[1:lenRest, ])
        names(newBlockGen) = names(blockGen(x))
        naFrameGen = as.data.frame(matrix(rep(NA, times = ncol(blockGen(x)) * nrow(newDf)), ncol = ncol(block(x))))
        names(naFrameGen) = names(oldBlockGen)
        restBlockGen = data.frame(oldBlockGen[-c(1:(lenRest + nrow(oldDf))), ])
        names(restBlockGen) = names(oldBlockGen)
        newBlockGen = rbind(newBlockGen, naFrameGen, restBlockGen)
        if (DB)
            print(newBlockGen)
        blockGen(x) = newBlockGen
    }
    oldBlock = block(x)
    newBlock = data.frame(oldBlock[1:lenRest, ])
    names(newBlock) = names(block(x))
    naFrame = as.data.frame(matrix(rep(max(block(x)[1:nrow(cube(x)), ]) + 1, times = ncol(block(x)) *
        nrow(newDf)), ncol = ncol(block(x))))
    names(naFrame) = names(oldBlock)
    restBlock = data.frame(oldBlock[-c(1:(lenRest + nrow(oldDf))), ])
    names(restBlock) = names(oldBlock)
    newBlock = rbind(newBlock, naFrame, restBlock)
    block(x) = newBlock
    x@centerStar <- newDf
    x
})
setGeneric("response", function(object) {
    standardGeneric("response")
})
setGeneric("response<-", function(object, value) {
    standardGeneric("response<-")
})
setMethod("response", "facDesign", function(object) {
    iIntern = order(runOrd(object))
    out = data.frame(object@response[iIntern, ])
    names(out) = names(object@response)
    return(out)
})
setReplaceMethod("response", "facDesign", function(object, value) {
    index = order(runOrd(object))
    if (!is.vector(value) && !is.data.frame(value))
        stop("vector or data.frame expected!")
    if (is.vector(value) && (is.numeric(value) || is.na(value))) {
        if (nrow(response(object)) != length(value))
            stop(paste("Number of rows for Design does not equal length of vector ", nrow(object),
                " != ", length(value), " "))
        object@response <- data.frame(value)
        object@response[index, ] <- value
        names(object@response) = make.names(deparse(substitute(value)))
        return(object)
    }
    if (is.data.frame(value)) {
        object@response <- value
        object@response[index, ] <- value
        object
    }
    return(object)
})
setGeneric("blockGen", function(object) {
    standardGeneric("blockGen")
})
setGeneric("blockGen<-", function(object, value) {
    standardGeneric("blockGen<-")
})
setMethod("blockGen", "facDesign", function(object) {
    return(object@blockGen)
})
setReplaceMethod("blockGen", "facDesign", function(object, value) {
    if (!is.vector(value) && !is.data.frame(value))
        stop("vector or data.frame expected!")
    if (is.vector(value) && (is.numeric(value) || is.na(value))) {
        if (nrow(object) != length(value))
            stop(paste("Number of rows for Design does not equal length of vector ", nrow(object),
                " != ", length(value), " "))
        object@blockGen <- as.data.frame(value)
        names(blockGen(object)) = deparse(substitute(value))
        object
    }
    if (is.data.frame(value)) {
        object@blockGen <- value
        object
    }
    return(object)
})
setGeneric("block", function(object) {
    standardGeneric("block")
})
setGeneric("block<-", function(object, value) {
    standardGeneric("block<-")
})
setMethod("block", "facDesign", function(object) {
    return(object@block)
})
setReplaceMethod("block", "facDesign", function(object, value) {
    if (!is.vector(value) && !is.data.frame(value))
        stop("vector or data.frame expected!")
    if (is.vector(value) && (is.numeric(value) || is.na(value))) {
        if (nrow(object) != length(value))
            stop(paste("Number of rows for Design does not equal length of vector ", nrow(object),
                " != ", length(value), " "))
        object@block <- as.data.frame(value)
        names(block(object)) = deparse(substitute(value))
        object
    }
    if (is.data.frame(value)) {
        object@block <- value
        object
    }
    return(object)
})
setGeneric("standOrd", function(x) standardGeneric("standOrd"))
setGeneric("standOrd<-", function(x, value) standardGeneric("standOrd<-"))
setMethod("standOrd", "facDesign", function(x) x@standardOrder)
setReplaceMethod("standOrd", "facDesign", function(x, value) {
    x@standardOrder <- value
    x
})
setGeneric("runOrd", function(x) standardGeneric("runOrd"))
setGeneric("runOrd<-", function(x, value) standardGeneric("runOrd<-"))
setMethod("runOrd", "facDesign", function(x) x@runOrder)
setReplaceMethod("runOrd", "facDesign", function(x, value) {
    x@runOrder <- value
    x
})
setGeneric(".generators", function(object) standardGeneric(".generators"))
setGeneric(".generators<-", function(x, value) standardGeneric(".generators<-"))
setMethod(".generators", "facDesign", function(object) object@generator)
setReplaceMethod(".generators", "facDesign", function(x, value) {
    x@generator <- value
    x
})
setGeneric("factors", function(x) standardGeneric("factors"))
setMethod("factors", "facDesign", function(x) x@factors)
setGeneric("factors<-", function(x, value) standardGeneric("factors<-"))
setReplaceMethod("factors", "facDesign", function(x, value) {
    if (length(value) != dim(cube(x))[2])
        stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
    x@factors <- value
    x
})
setMethod("names", "facDesign", function(x) {
    return(sapply(x@factors, names))
})
setReplaceMethod("names", "facDesign", function(x, value) {
    for (i in 1:length(x@factors)) names(x@factors[[i]]) = as.character(value[i])
    x
})
setMethod("as.data.frame", "facDesign", function(x, row.names = NULL, optional = FALSE,
    ...) {
    if (!is.null(cube(x))) {
        frameOut = cube(x)
    }
    else return(NULL)
    if (!is.null(centerCube))
        frameOut = rbind(frameOut, centerCube(x))
    if (!is.null(star(x)))
        frameOut = rbind(frameOut, star(x))
    if (!is.null(centerStar(x)))
        frameOut = rbind(frameOut, centerStar(x))
    if (!is.null(factors(x)) && length(factors(x)) == dim(frameOut)[2]) {
        names(frameOut) = as.character(names(names(x)))
    }
    if (nrow(blockGen(x)) > 0) {
    }
    if (nrow(block(x)) > 0) {
        frameOut = cbind(block(x), frameOut)
    }
    if (nrow(runOrd(x)) > 0) {
        frameOut = cbind(runOrd(x), frameOut)
    }
    if (nrow(standOrd(x)) > 0) {
        frameOut = cbind(standOrd(x), frameOut)
    }
    if (nrow(frameOut) == nrow(response(x)))
        frameOut = cbind(frameOut, x@response)
    else {
        temp = as.data.frame(matrix(NA, nrow = nrow(frameOut), ncol = ncol(response(x))))
        names(temp) = names(response(x))
        frameOut = cbind(frameOut, temp)
    }
    runIndex = order(runOrd(x))
    out = frameOut[runIndex, ]
    return(out)
})
as.data.frame.facDesign = function(x, row.names = NULL, optional = FALSE, ...) {
    if (!is.null(cube(x))) {
        frameOut = cube(x)
    }
    else return(NULL)
    if (!is.null(centerCube))
        frameOut = rbind(frameOut, centerCube(x))
    if (!is.null(star(x)))
        frameOut = rbind(frameOut, star(x))
    if (!is.null(centerStar(x)))
        frameOut = rbind(frameOut, centerStar(x))
    if (!is.null(factors(x)) && length(factors(x)) == dim(frameOut)[2]) {
        names(frameOut) = as.character(names(names(x)))
    }
    if (nrow(blockGen(x)) > 0) {
    }
    if (nrow(block(x)) > 0) {
        frameOut = cbind(block(x), frameOut)
    }
    if (nrow(runOrd(x)) > 0) {
        frameOut = cbind(runOrd(x), frameOut)
    }
    if (nrow(standOrd(x)) > 0) {
        frameOut = cbind(standOrd(x), frameOut)
    }
    if (nrow(frameOut) == nrow(response(x)))
        frameOut = cbind(frameOut, x@response)
    else {
        temp = as.data.frame(matrix(NA, nrow = nrow(frameOut), ncol = ncol(response(x))))
        names(temp) = names(response(x))
        frameOut = cbind(frameOut, temp)
    }
    runIndex = order(runOrd(x))
    out = frameOut[runIndex, ]
    return(frameOut)
}
setMethod("show", signature(object = "facDesign"), function(object) {
    runIndex = order(runOrd(object))
    print(format(as.data.frame(object), digits = 4))
    invisible(as.data.frame(object))
})
setMethod("summary", signature(object = "facDesign"), function(object) {
    doeFactors = factors(object)
    cat("Information about the factors:\n\n")
    .nfp(object)
    cat("-----------\n")
    print(as.data.frame(object))
    temp = aliasTable(object, show = FALSE)
    if (ncol(temp) > 0) {
        cat("\n---------\n\n")
        identity(object)
        cat("\n")
    }
    invisible(as.data.frame(object))
})
.identityOld = function(x, DB = FALSE) {
    identityList = vector(mode = "list", length = 0)
    varName = deparse(substitute(x))
    ident = numeric(0)
    resolution = numeric(0)
    x = cube(x)
    index = numeric(0)
    for (i in 1:(dim(x)[2])) {
        if (!(TRUE && (unique(x[, i]) %in% c(-1, 1))))
            index = c(index, i)
    }
    if (length(index) > 0) {
        x = x[, -index]
        cat("Column(s) ", index, " are discarded for analysis\n")
    }
    n = dim(x)[2]
    if (n <= 1)
        stop("Factorial Design contains only one row!")
    a = 0
    for (m in 2:n) {
        combMat = combn(1:n, m)
        for (i in 1:(dim(combMat)[2])) {
            ident = NULL
            temp = x[, combMat[, i]]
            colProd = apply(temp, 1, prod)
            if (length(unique(colProd)) == 1) {
                if (DB) {
                  cat("\n")
                  cat("Identitaet gefunden\n")
                }
                if (sum(colProd) < 0)
                  ident = -combMat[, i]
                else ident = combMat[, i]
                a = a + 1
                identityList[[a]] = ident
                charIdentity = character(0)
                for (j in 1:length(ident)) charIdentity = c(charIdentity, names(x)[ident[j]])
                names(identityList)[[a]] = paste(charIdentity, sep = "", collapse = "")
            }
            if (DB) {
                cat("ident: ", ident, "\n")
                print(combMat[, i])
                print(apply(temp, 1, prod))
            }
        }
    }
    cat("Defining relations:\n")
    if (length(identityList) > 0) {
        for (i in 1:length(identityList)) {
            identLen = length((strsplit(names(identityList)[i], split = character(0))[[1]]))
            if (length(resolution) == 0 || identLen > resolution)
                resolution = c(resolution, identLen)
            cat("I = ", names(identityList)[i], "\t\tColumns:", identityList[[i]], "\n")
        }
        cat("\nResolution: ", as.character(as.roman(min(resolution))), "\n")
    }
    invisible(identityList)
}
.helpAliasTable = function(fdo, k) {
    X = unique(cube(fdo))
    N = nrow(X)
    columns = names(X[, 1:k])
    X1 = matrix(1, nrow = N, ncol = 1)
    nameVec = c("Identity")
    for (i in 1:k) {
        temp = combn(columns, i)
        for (j in 1:ncol(temp)) {
            index = names(names(fdo)) %in% temp[, j]
            if (length((1:length(index))[index]) == 1) {
                X1 = cbind(X1, X[, index])
                nameVec = c(nameVec, temp[, j])
            }
            else {
                X1 = cbind(X1, apply(X[, index], 1, prod))
                nameVec = c(nameVec, paste(temp[, j], sep = "", collapse = ""))
            }
        }
        X1 = data.frame(X1)
        names(X1) = nameVec
    }
    return(X1)
}
setGeneric("identity")
setMethod("identity", signature(x = "facDesign"), function(x) {
    identity = character(0)
    identityList = vector(mode = "list", length = 0)
    resolution = numeric(0)
    temp = NULL
    A = aliasTable(x, show = FALSE)
    if (any(dim(A) == 0))
        return(identityList)
    temp = as.matrix(A["Identity", ])
    boolTemp = apply(temp, 2, as.logical)
    identity = row.names(temp)[boolTemp[, 1]]
    if (length(identity) > 0) {
        charList = strsplit(toupper(identity), split = "")
        identityList = lapply(charList, match, .NAMES[1:25])
        names(identityList) = identity
    }
    cat("Defining relations:\n")
    if (length(identityList) > 0) {
        for (i in 1:length(identityList)) {
            identLen = length((strsplit(names(identityList)[i], split = character(0))[[1]]))
            if (length(resolution) == 0 || identLen > resolution)
                resolution = c(resolution, identLen)
            cat("I = ", names(identityList)[i], "\t\tColumns:", identityList[[i]], "\n")
        }
        cat("\nResolution: ", as.character(as.roman(min(resolution))), "\n")
    }
    invisible(identityList)
})
confounds = function(x, depth = 2) {
    DB = FALSE
    varName = deparse(substitute(x))
    identityList = identity(x)
    x = cube(x)
    if (length(identityList) < 1) {
        print(paste(varName, " contains no defining relations!"))
        invisible()
    }
    effect1 = numeric(0)
    effect2 = numeric(0)
    if (DB)
        identityList
    index = numeric(0)
    for (i in 1:(dim(x)[2])) {
        if (!(TRUE && (unique(x[, i]) %in% c(-1, 1))))
            index = c(index, i)
    }
    if (length(index) > 0)
        x = x[, -index]
    if (DB)
        cat("Column(s) ", index, " are discarded for analysis\n")
    n = dim(x)[2]
    if (n <= 1)
        stop("Factorial Design contains only one row!")
    for (j in 1:length(identityList)) {
        ident = identityList[[j]]
        for (m in 1:n) {
            combMat = combn(1:n, m)
            for (i in 1:(dim(combMat)[2])) {
                isect = intersect(ident, combMat[, i])
                conf = setdiff(ident, isect)
                conf = sort(c(conf, setdiff(combMat[, i], isect)))
                effect1 = c(effect1, paste(sort(names(x)[as.numeric((combMat[, i]))]), sep = "",
                  collapse = ""))
                effect2 = c(effect2, paste(sort(names(x)[conf]), sep = "", collapse = ""))
                if (DB) {
                  cat("Effect(s) ", as.numeric((combMat[, i])), " aliased with Effect(s)", conf,
                    "\n")
                  cat("Effect(s)", names(x)[as.numeric((combMat[, i]))], " aliased with Effects ",
                    names(x)[conf], "\n")
                }
            }
        }
    }
    if (DB) {
        cat(effect1, "\n")
        cat(effect2, "\n")
    }
    if (length(effect1) > 0)
        dupIndex = numeric(0)
    for (i in 1:length(effect1)) {
        if (DB) {
            cat("i: ", i, "\tlength(effect1): ", length(effect1), "\n")
            cat("effect 1 : ", effect1, "\n")
        }
        if (i > length(effect1))
            break
        index = (1:length(effect1))[effect2 == effect1[i]]
        if (DB)
            cat("index: ", index, "\n")
        dupIndex = numeric(0)
        for (j in index) {
            if (effect1[j] == effect2[i]) {
                if (i != j)
                  dupIndex = c(dupIndex, j)
            }
        }
        if (length(dupIndex > 0)) {
            effect1 = effect1[-dupIndex]
            effect2 = effect2[-dupIndex]
        }
    }
    cat("\nAlias Structure:\n")
    for (i in 1:length(effect1)) {
        if ((length(strsplit(effect1[i], split = character(0))[[1]]) <= depth) && (length(strsplit(effect2[i],
            split = character(0))[[1]]) <= depth))
            cat(effect1[i], "\tis confounded with\t", effect2[i], "\n")
        if (identical(depth, "all"))
            cat(effect1[i], "\tis confounded with\t", effect2[i], "\n")
    }
    invisible(effect1)
}
fracDesign = function(k = 3, p = 0, gen = NULL, replicates = 1, blocks = 1, centerCube = 0, random.seed = 1234) {    ###
    DB = FALSE
    STDfdo=FALSE
    if (p<0 || p>7)                                                                                                  ###
     stop("p needs to be an integer between 0 and 7!")                                                               ###
    if (abs(p - round(p)) > .Machine$double.eps^0.5)                                                                 ###
    {                                                                                                                ###
     warning(paste("p needs to be an integer but is real! p was rounded to", round(p)))                              ###
     p=round(p)                                                                                                      ###
    }                                                                                                                ###
    if (p != 0)                                                                                                      ###
    {                                                                                                                ###
     gen = NULL                                                                                                      ###
     for (i in 1:length(.fdoOrth))                                                                                   ###
      {                                                                                                              ###
       if (k==.fdoOrth[[i]]$k && p==.fdoOrth[[i]]$p)                                                                 ###
       {                                                                                                             ###
         STDfdo=TRUE                                                                                                 ###
         return(fracDesign(k=.fdoOrth[[i]]$k, gen=.fdoOrth[[i]]$gen,                                                 ###
                replicates = replicates, blocks = blocks,                                                            ###
                centerCube = centerCube, random.seed = random.seed))                                                 ###
       }                                                                                                             ###
      }                                                                                                              ###
      if(STDfdo==FALSE)                                                                                              ###
       stop("No standard Design for the choosen combination of k and p (see: fracChoose())!")                        ###
    }                                                                                                                ###
    if (!is.numeric(random.seed))
        stop("random.seed needs to be numeric")
    if (!is.numeric(blocks))
        stop("blocks needs to be numeric!")
    if (!is.numeric(replicates))
        stop("replicates needs to be numeric!")
    else if (replicates < 0)
        stop("replicates need to >= 0")
    N <- 2^k
    X <- matrix(NA, nrow = N, ncol = k)
    for (j in 1:k) X[, j] <- rep(sort(rep(c(-1, 1), N/2^j)), 2^(j - 1))
    X <- X[, ncol(X):1]
    if (is.null(gen)) {
        X = as.data.frame(X)
        names(X) = .NAMES[1:k]
    }
    origX = X
    if (replicates > 1) {
        for (i in 2:replicates) {
            X = rbind(X, origX)
        }
    }
    frameOut = data.frame(X)
    if (DB)
        print("juhu")
    if (!is.null(gen)) {
        listGen = vector("list", length(gen))
        .numFactors = numeric(0)
        charFactors = character(0)
        if (DB) {
            cat(paste("gen is: ", gen, "\n"))
            cat(paste("length of gen is: ", length(gen), "\n"))
            print(listGen)
        }
        temp = character(0)
        for (i in seq(along = gen)) {
            if (DB)
                cat("gen[", i, "] = ", gen[i], "\n")
            if (!is.character(gen[i]))
                stop("Defining Relations should contain characters only!")
            chars = strsplit(gen[i], split = character(0))[[1]]
            if (DB) {
                cat("\nchars: ")
                print(chars)
                cat("\n")
            }
            checkDupl = character(0)
            for (j in 1:length(chars)) {
                if (chars[j] %in% toupper(c(.NAMES[1:26], letters[1:26]))) {
                  if (chars[j] %in% checkDupl)
                    stop("Defining relations contain one or more duplicates!")
                  checkDupl = c(checkDupl, chars[j])
                  temp = c(temp, chars[j])
                }
            }
            if (DB) {
                cat("\ntemp: ")
                print(temp)
                cat("\n")
            }
        }
        temp = sort(unique(temp))
        numCharVec = 1:length(temp)
        names(numCharVec) = temp
        if (DB) {
            cat("Zuordnung Buchstabe und Spalte:\n")
            print(numCharVec)
            cat("\n")
        }
        for (i in seq(along = gen)) {
            if (DB)
                cat("gen[", i, "] = ", gen[i], "\n")
            if (!is.character(gen[i]))
                stop("Defining Relations should contain characters only!")
            chars = strsplit(gen[i], split = character(0))[[1]]
            numVec = numeric(0)
            charVec = character(0)
            allowedChars = c(.NAMES[1:26], letters[1:26], "=")
            for (j in 1:length(chars)) {
                if (chars[j] %in% allowedChars) {
                  if ((chars[j] == "=") & (length(numVec) != 1))
                    stop("check position of \"=\" in generators!")
                  if (chars[j] != "=") {
                    charVec = c(charVec, toupper(chars[j]))
                    numVec = c(numVec, numCharVec[names(numCharVec) == toupper(chars[j])])
                  }
                }
            }
            if (DB) {
                cat("charVec for i = ", i, ": ", charVec, "\n")
                cat("numVec for i = ", i, ": ", numVec, "\n")
            }
            listGen[[i]] = numVec
            .numFactors = c(.numFactors, numVec)
            charFactors = c(charFactors, charVec)
        }
        if (DB)
            print("juhu")
        names(.numFactors) = charFactors
        if (length(unique(.numFactors)) > k)
            stop("number of distinct Factors in generators greater than k!")
        if (DB) {
            print(listGen)
            print(.numFactors)
            print(charFactors)
        }
        for (i in seq(along = listGen)) {
            ind <- trunc(listGen[[i]])
            if (any(abs(ind) > k))
                stop(paste("generator:", paste(ind[1], "=", paste(ind[-1], collapse = "*")), "includes undefined columns"))
            x <- rep(sign(ind[1]), N)
            for (j in ind[-1]) x <- x * X[, j]
            X[, abs(ind[1])] <- x
        }
        X <- unique(X)
        origX = X
        if (replicates > 1) {
            for (i in 2:replicates) {
                X = rbind(X, origX)
            }
        }
        frameOut = as.data.frame(X)
        names(frameOut) = names(numCharVec)
        if (k > length(temp)) {
            charsLeft = (.NAMES[1:26])[-match(charFactors, .NAMES[1:26])]
            naIndex = (1:k)[is.na(names(frameOut))]
            names(frameOut)[naIndex] = charsLeft[1:length(naIndex)]
        }
    }
    DesignOut = new("facDesign")
    DesignOut@generator = gen
    cube(DesignOut) = frameOut
    listFac = vector("list", ncol(frameOut))
    for (i in seq(along = listFac)) listFac[i] = new("doeFactor")
    names(listFac) = names(frameOut)
    factors(DesignOut) = listFac
    if (DB)
        print(frameOut)
    if (DB)
        print("yes")
    if (DB)
        print("aha")
    numRows = nrow(cube(DesignOut)) + nrow(star(DesignOut)) + nrow(centerStar(DesignOut)) + nrow(centerCube(DesignOut))
    if (DB) {
        print(numRows)
        print("response")
    }
    DesignOut@response = data.frame(y = rep(NA, numRows))
    if (DB)
        print("response")
    standardOrder = data.frame(matrix(data = 1:numRows, nrow = numRows, ncol = 1))
    names(standardOrder) = "StandOrder"
    standardOrder
    standOrd(DesignOut) = standardOrder
    if (DB)
        print("1")
    set.seed(random.seed)
    runOrder = as.data.frame(standardOrder[sample(1:numRows), ])
    if (DB)
        print("2")
    names(runOrder) = "RunOrder"
    runOrd(DesignOut) = runOrder
    if (DB)
        print("3")
    if (centerCube >= 1) {
        temp = data.frame(matrix(rep(0, centerCube * k), ncol = k, nrow = centerCube))
        names(temp) = names(frameOut)
        centerCube(DesignOut) = temp
    }
    temp = try(blocking(DesignOut, blocks = blocks))
    if (inherits(temp, "try-error"))
        stop("Blocking not possible!")      #return(DesignOut)                  ###
    return(blocking(DesignOut, blocks = blocks))
}
facDesign = function(k = 3, p = 0, replicates = 1, blocks = 1, centerCube = 0) {  ###
    frameOut = fracDesign(k = k, p = p, gen = NULL, replicates = replicates, blocks = blocks, centerCube = centerCube)  ###
    return(frameOut)
}

###

.fdoOrth = vector(mode = "list", length = 3)                                                           ###
.fdoOrth[[1]] = list(k = 3, gen = "C=AB", p = 1)                                                       ###
.fdoOrth[[2]] = list(k = 4, gen = "D=ABC", p = 1)                                                      ###
.fdoOrth[[3]] = list(k = 5, gen = c("D=AB","E=AC"), p = 2)                                             ###
.fdoOrth[[4]] = list(k = 6, gen = c("D=AB","E=AC","F=BC"), p = 3)                                      ###
.fdoOrth[[5]] = list(k = 7, gen = c("D=AB","E=AC","F=BC","G=ABC"), p = 4)                              ###
.fdoOrth[[6]] = list(k = 5, gen = "E=ABCD", p = 1)                                                     ###
.fdoOrth[[7]] = list(k = 6, gen = c("E=ABC","F=BCD"), p = 2)                                           ###
.fdoOrth[[8]] = list(k = 7, gen = c("E=ABC","F=BCD","G=ACD"), p = 3)                                   ###
.fdoOrth[[9]] = list(k = 8, gen = c("E=BCD","F=ACD","G=ABC","H=ABD"), p = 4)                           ###
.fdoOrth[[10]] = list(k = 9, gen = c("E=ABC","F=BCD","G=ACD","H=ABD","J=ABCD"), p = 5)                 ###
.fdoOrth[[11]] = list(k = 10, gen = c("E=ABC","F=BCD","G=ACD","H=ABD","J=ABCD","K=AB"), p = 6)         ###
.fdoOrth[[12]] = list(k = 11, gen = c("E=ABC","F=BCD","G=ACD","H=ABD","J=ABCD","K=AB","L=AC"), p = 7)  ###
.fdoOrth[[13]] = list(k = 6, gen = "F=ABCDE", p = 1)                                                   ###
.fdoOrth[[14]] = list(k = 7, gen = c("F=ABCD","G=ABDE"), p = 2)                                        ###
.fdoOrth[[15]] = list(k = 8, gen = c("F=ABC","G=ABD","H=BCDE"), p = 3)                                 ###
.fdoOrth[[16]] = list(k = 9, gen = c("F=BCDE","G=ACDE","H=ABDE","J=ABCE"), p = 4)                      ###
.fdoOrth[[17]] = list(k = 10, gen = c("F=ABCD","G=ABCE","H=ABDE","J=ACDE","K=BCDE"), p = 5)            ###
.fdoOrth[[18]] = list(k = 11, gen = c("F=ABC","G=BCD","H=CDE","J=ACD","K=AEF","L=ADEF"), p = 6)        ###
.fdoOrth[[19]] = list(k = 7, gen = "G=ABCDEF", p = 1)                                                  ###
.fdoOrth[[20]] = list(k = 8, gen = c("G=ABCD","H=ABEF"), p = 2)                                        ###
.fdoOrth[[21]] = list(k = 9, gen = c("G=ABCD","H=ABEF","J=CDEF"), p = 3)                               ###
.fdoOrth[[22]] = list(k = 10, gen = c("G=BCDF","H=ACDF","J=ABDE","K=ABCE"), p = 4)                     ###
.fdoOrth[[23]] = list(k = 11, gen = c("G=CDE","H=ABCD","J=ABF","K=BDEF","L=ADEF"), p = 5)              ###
.fdoOrth[[24]] = list(k = 8, gen = "H=ABCDEFG", p = 1)                                                 ###
.fdoOrth[[25]] = list(k = 9, gen = c("H=ACDFG","J=BCEFG"), p = 2)                                      ###
.fdoOrth[[26]] = list(k = 10, gen = c("H=ABCG","J=BCDE","K=ACDF"), p = 3)                              ###
.fdoOrth[[27]] = list(k = 11, gen = c("H=ABCG","J=BCDE","K=ACDF","L=ABCDEFG"), p = 4)                  ###