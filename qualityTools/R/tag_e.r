.NAMES = LETTERS[c(1:8, 10:26)]
setClass(Class = "taguchiDesign", representation = representation(name = "character", factors = "list", design = "data.frame", designType = "character", 
    replic = "data.frame", response = "data.frame", Type = "data.frame", block = "data.frame", runOrder = "data.frame", standardOrder = "data.frame", desireVal = "list", 
    desirability = "list", fits = "data.frame"))
setClass("taguchiFactor", representation = representation(values = "ANY", name = "character", unit = "character", type = "character"), prototype = prototype(values = NA, 
    name = " ", unit = " ", type = "numeric"))
setGeneric("values", function(object) standardGeneric("values"))
setGeneric("values<-", function(object, value) standardGeneric("values<-"))
setMethod("values", "taguchiDesign", function(object) {
    listOut = vector(mode = "list")
    for (i in names(object@design)) {
        listOut[[i]] = .values(object@factors[[i]])
    }
    return(listOut)
})
setReplaceMethod("values", "taguchiDesign", function(object, value) {
    for (i in names(value)) {
        if (i %in% names(object@design)) 
            if (length(value[[i]]) == length(unique(object@design[, i]))) 
                .values(object@factors[[i]]) = value[[i]]
            else stop("Number of values greater or less than number of factor settings!")
    }
    return(object)
})
setGeneric(".values", function(object) standardGeneric(".values"))
setGeneric(".values<-", function(object, value) standardGeneric(".values<-"))
setMethod(".values", "taguchiFactor", function(object) object@values)
setReplaceMethod(".values", "taguchiFactor", function(object, value) {
    object@values <- value
    return(object)
})
setMethod(".unit", "taguchiFactor", function(object) object@unit)
setReplaceMethod(".unit", "taguchiFactor", function(x, value) {
    x@unit <- value
    x
})
setMethod("units", "taguchiDesign", function(x) {
    return(sapply(factors(x), .unit))
})
setMethod("units<-", "taguchiDesign", function(x, value) {
    for (i in 1:length(x@factors)) if (length(value) > 1) 
        .unit(x@factors[[i]]) = as.character(value[i])
    else .unit(x@factors[[i]]) = as.character(value[1])
    x
})
setMethod("factors", "taguchiDesign", function(x) x@factors)
setReplaceMethod("factors", "taguchiDesign", function(x, value) {
    if (length(value) != ncol(x@design)) 
        stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
    x@factors <- value
    x
})
setMethod("names", "taguchiDesign", function(x) {
    return(sapply(x@factors, names))
})
setReplaceMethod("names", "taguchiDesign", function(x, value) {
    for (i in 1:length(x@factors)) names(x@factors[[i]]) = as.character(value[i])
    x
})
setMethod("names", "taguchiFactor", function(x) {
    x@name
})
setReplaceMethod("names", "taguchiFactor", function(x, value) {
    x@name <- value
    x
})
setMethod("as.data.frame", "taguchiDesign", function(x, row.names = NULL, optional = FALSE) {
    frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
    return(frameOut)
})
as.data.frame.taguchiDesign = function(x, row.names = NULL, optional = FALSE, ...) {
    frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
    return(frameOut)
}
setMethod("show", signature(object = "taguchiDesign"), function(object) {
    print(format(as.data.frame(object), digits = 4))
})
setMethod("response", "taguchiDesign", function(object) {
    return(object@response)
})
setReplaceMethod("response", "taguchiDesign", function(object, value) {
#    print(deparse(substitute(value)))                                          ###
    if (!is.numeric(value) & !is.data.frame(value)) 
        stop("vector or data.frame must be given")
    if (is.numeric(value)) {
        if (length(value) != nrow(object@design)) 
            stop("differing lengths")
        temp = data.frame(value)
        names(temp) = deparse(substitute(value))[1]
        value = temp
    }
    if (is.data.frame(value)) {
        if (nrow(value) != nrow(object@design)) 
            stop("differing number of rows")
    }
    object@response = value
    return(object)
})
setMethod(".nfp", "taguchiDesign", function(object) {
    x = factors(object)
    DB = FALSE
    if (is.list(x) && length(x[[1]]) > 0) {
        numAttr = length(attributes(x[[1]])) - 1
        .numFac = length(x)
        len = 0
        for (i in names(x)) if (length(x[[i]]@values) > len) 
            len = length(x[[i]]@values)
        numAttr = numAttr + len
        numrows = numAttr - 1
        frameOut = data.frame(matrix(NA, ncol = .numFac, nrow = numrows))
        names(frameOut) = names(x)
        rownames(frameOut) = c(paste("value", 1:len), "name", "unit", "type")
        for (i in names(x)) {
            vin = 1:length(x[[i]]@values)
            frameOut[vin, i] = x[[i]]@values
            frameOut[numrows - 2, i] = x[[i]]@name
            frameOut[numrows - 1, i] = x[[i]]@unit
            frameOut[numrows, i] = x[[i]]@type
        }
        print(frameOut)
    }
})
setMethod("show", signature(object = "taguchiDesign"), function(object) {
    print(as.data.frame(object))
})
setMethod("summary", signature(object = "taguchiDesign"), function(object) {
    cat(paste("Taguchi", toupper(object@designType), "Design"))
    cat("\n")
    cat("Information about the factors:\n\n")
    .nfp(object)
    cat("\n")
    cat("-----------\n")
    cat("\n")
    print(as.data.frame(object))
    cat("\n")
    cat("-----------\n")
    cat("\n")
})
taguchiDesign = function(design, randomize = TRUE, replicates = 1) {
    DB = FALSE
    odo = NA
    type = "single"
    for (i in seq(along = .oaList)) {
        pmatch(design, .oaList[[i]]$id)
        if (!is.na(pmatch(design, .oaList[[i]]$id))) {
            if (DB) 
                print(design)
            temp = .oaList[[i]]
            design = temp$design
            repVec = rep(1, nrow(design))
            if (replicates > 1) {
                X = temp$design
                for (i in 1:(replicates - 1)) {
                  design = rbind(design, X)
                  repVec = c(repVec, rep(i + 1, times = nrow(X)))
                }
            }
            Replicate = data.frame(Replicate = as.numeric(repVec))
            if (DB) 
                print(Replicate)
            odo = new("taguchiDesign")
            odo@design = design
            names(odo@design) = .NAMES[1:ncol(design)]
            odo@name = temp$id
            odo@designType = temp$type
            odo@replic = Replicate
            StandOrder = 1:nrow(odo@design)
            RunOrder = StandOrder
            if (randomize) {
                RunOrder = sample(1:nrow(odo@design), nrow(odo@design), replace = FALSE, prob = NULL)
            }
            odo@design = odo@design[order(RunOrder), ]
            odo@replic = data.frame(Replicate = odo@replic[order(RunOrder), 1])
            row.names(odo@design) = odo@design$RunOrder
            odo@runOrder = data.frame(RunOrder = data.frame(RunOrder = RunOrder)[order(RunOrder), ])
            odo@standardOrder = data.frame(StandOrder = data.frame(StandOrder = StandOrder)[order(RunOrder), ])
            odo@response = data.frame(y = rep(NA, nrow(odo@design)))
            tfList = vector("list", ncol(design))
            for (i in seq(along = tfList)) tfList[i] = new("taguchiFactor")
            names(tfList) = names(odo@design)
            factors(odo) = tfList
            valList = list(length = length(names(odo)))
            for (i in names(names(odo))) valList[[i]] = sort(unique(odo@design[, i]))
            values(odo) = valList
            return(odo)
        }
    }
    return(NA)
}
oaChoose = function(factors1, factors2, level1, level2, ia) {
    params = list(factors1 = 0, factors2 = 0, level1 = 0, level2 = 0, ia = 0)
    if (!missing(ia)) 
        params$ia = ia
    if (!missing(factors2)) 
        params$factors2 = factors2
    if (!missing(level2)) 
        params$level2 = level2
    do.call(taguchiChoose, params)
}
taguchiChoose = function(factors1 = 0, factors2 = 0, level1 = 0, level2 = 0, ia = 0) {
    if (factors1 == 0 & factors2 == 0 & level1 == 0 & level2 == 0 & ia == 0) {
        temp = vector(mode = "character", length = length(.oaList))
        for (i in 1:length(.oaList)) temp[i] = .oaList[[i]]$id
        temp = c(temp, rep(" ", (length(temp)%/%6 + 1) * 6 - length(temp)))
        mat = data.frame(matrix(temp, ncol = 6, byrow = TRUE))
        names(mat) = rep(" ", ncol(mat))
        print(mat)
        cat("\n")
        cat("Choose a design using e.g. taguchiDesign(\"L4_2\")")
        cat("\n")
    }
    else {
        DB = FALSE
        if (factors2 <= 0) 
            level2 = 0
        if (DB) 
            print(level2)
        Anzahl_Spalten = factors1 + factors2 + ia
        ss = list()
        for (i in seq(along = .oaList)) {
            li = .oaList[[i]]
            if (li$factors1 >= factors1 & li$factors2 >= factors2 & (li$levels1 == level1 | li$levels1 == level2) & (li$levels2 == level2 | li$levels2 == level1) & 
                li$anzahl_spalten >= Anzahl_Spalten) 
                ss[i] = li$id
        }
        out = as.character(ss)
        out = out[out != "NULL"]
        if (length(out) > 0) {
            cat(paste(factors1, "factors on", level1, "levels and", factors2, "factors on", level2, "levels with", ia, "desired interactions to be estimated\n"))
            cat("\n")
            cat("Possible Designs:\n")
            cat("\n")
            cat(paste(out, sep = " | "))
            cat("\n")
            cat("\n")
            cat(paste("Use taguchiDesign(\"", out[1], "\") or different to create a taguchi design object\n", sep = ""))
        }
        else {
            cat("No Design Found\n")
            cat("\n")
            out = NA
        }
        invisible(out)
    }
}
.replace2s = function(x) {
    if (!is.data.frame(x)) 
        stop(paste(deparse(substitute(x)), "needs to be a data.frame"))
    for (i in 1:ncol(x)) x[x[, i] == 2, i] = -1
    return(x)
}
.helpAliasTable = function(fdo, k, degree = 3) {
    if (degree > k) {
        degree = k
    }
    if (class(fdo) == "facDesign") 
        X = unique(cube(fdo))
    if (class(fdo) == "taguchiDesign") {
        X = unique(fdo@design)
        X = .replace2s(X)
    }
    N = nrow(X)
    columns = names(X[, 1:k])
    X1 = matrix(1, nrow = N, ncol = 1)
    nameVec = c("Identity")
    for (i in 1:degree) {
        temp = combn(columns, i)
        for (j in 1:ncol(temp)) {
            if (class(fdo) == "facDesign") 
                index = names(names(fdo)) %in% temp[, j]
            if (class(fdo) == "taguchiDesign") 
                index = names(X) %in% temp[, j]
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
aliasTable = function(fdo, degree, show = TRUE) {
    if (class(fdo) == "facDesign") {
        X = unique(cube(fdo))
        N = nrow(X)
        k = log2(N)
        kPlusP = ncol(X)
        if (missing(degree)) 
            degree = min(c(4, k + 1))
        X1 = .helpAliasTable(fdo, k, degree = degree - 1)
        X2 = .helpAliasTable(fdo, k = kPlusP, degree)
    }
    if (class(fdo) == "taguchiDesign") {
        if (length(table(as.numeric(as.matrix(fdo@design)))) != 2) 
            stop("calculation of an alias table for mixed designs is not supported")
        k = ncol(fdo@design)
        if (missing(degree)) 
            degree = min(c(3, k))
        X1 = unique(fdo@design)
        X1 = .replace2s(X1)
        X2 = .helpAliasTable(fdo, k, degree)
        X1 = cbind(data.frame(Identity = rep(1, times = nrow(X1))), X1)
    }
    logVec = !(names(X2) %in% names(X1))
    X2 = X2[, logVec]
    X1 = as.matrix(X1)
    X2 = as.matrix(X2)
    alias.matrix = solve(t(X1) %*% X1) %*% t(X1) %*% X2
    if (show) 
        print(round(alias.matrix, 2))
    invisible(alias.matrix)
}
setMethod("identity", signature(x = "taguchiDesign"), function(x) {
    identity = character(0)
    identityList = vector(mode = "list", length = 0)
    resolution = numeric(0)
    temp = NULL
    A = aliasTable(x)
    if (any(dim(A) == 0)) 
        return(identityList)
    temp = as.matrix(A["Identity", ])
    boolTemp = apply(temp, 2, as.logical)
    identity = row.names(temp)[boolTemp[, 1]]
    if (length(identity) > 0) {
        charList = strsplit(toupper(identity), split = "")
        identityList = lapply(charList, match, LETTERS[1:26])
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

################################################################################
#.NAMES = LETTERS[c(1:8, 10:26)]
setClass(Class = "pbDesign", representation = representation(name = "character", factors = "list", design = "data.frame", designType = "character",
    replic = "data.frame", response = "data.frame", Type = "data.frame", block = "data.frame", runOrder = "data.frame", standardOrder = "data.frame", desireVal = "list",
    desirability = "list", fits = "data.frame"))

setClass("pbFactor", representation = representation(values = "ANY", name = "character", unit = "character", type = "character"), prototype = prototype(values = NA,
    name = " ", unit = " ", type = "numeric"))

#setGeneric("values", function(object) standardGeneric("values"))
#setGeneric("values<-", function(object, value) standardGeneric("values<-"))
setMethod("values", "pbDesign", function(object) {
    listOut = vector(mode = "list")
    for (i in names(object@design)) {
        listOut[[i]] = .values(object@factors[[i]])
    }
   return(listOut)
})
setReplaceMethod("values", "pbDesign", function(object, value) {
    for (i in names(value)) {
        if (i %in% names(object@design))
            if (length(value[[i]]) == length(unique(object@design[, i])))
                .values(object@factors[[i]]) = value[[i]]
            else stop("Number of values greater or less than number of factor settings!")
    }
    return(object)
})

#setGeneric(".values", function(object) standardGeneric(".values"))
#setGeneric(".values<-", function(object, value) standardGeneric(".values<-"))
setMethod(".values", "pbFactor", function(object) object@values)
setReplaceMethod(".values", "pbFactor", function(object, value) {
    object@values <- value
    return(object)
})
setMethod(".unit", "pbFactor", function(object) object@unit)
setReplaceMethod(".unit", "pbFactor", function(x, value) {
    x@unit <- value
    x
})
setMethod("names", "pbFactor", function(x) {
    x@name
})
setReplaceMethod("names", "pbFactor", function(x, value) {
    x@name <- value
    x
})

setMethod("units", "pbDesign", function(x) {
    return(sapply(factors(x), .unit))
})
setMethod("units<-", "pbDesign", function(x, value) {
    for (i in 1:length(x@factors)) if (length(value) > 1)
        .unit(x@factors[[i]]) = as.character(value[i])
    else .unit(x@factors[[i]]) = as.character(value[1])
    x
})
setMethod("factors", "pbDesign", function(x) x@factors)
setReplaceMethod("factors", "pbDesign", function(x, value) {
    if (length(value) != ncol(x@design))
        stop("\nNumber of factors doesn't match with number of columns for factorial Design\n")
    x@factors <- value
    x
})
setMethod("names", "pbDesign", function(x) {
    return(sapply(x@factors, names))
})
setReplaceMethod("names", "pbDesign", function(x, value) {
    for (i in 1:length(x@factors)) names(x@factors[[i]]) = as.character(value[i])
    x
})
setMethod("names", "pbFactor", function(x) {
    x@name
})

setMethod("as.data.frame", "pbDesign", function(x, row.names = NULL, optional = FALSE) {
    frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
    return(frameOut)
})
as.data.frame.pbDesign = function(x, row.names = NULL, optional = FALSE, ...) {
    frameOut = cbind(x@standardOrder, x@runOrder, x@replic, x@design, x@response)
    return(frameOut)
}
setMethod("show", signature(object = "pbDesign"), function(object) {
    print(format(as.data.frame(object), digits = 4))
})
setMethod("response", "pbDesign", function(object) {
    return(object@response)
})
setReplaceMethod("response", "pbDesign", function(object, value) {
#    print(deparse(substitute(value)))                                          ###
    if (!is.numeric(value) & !is.data.frame(value))
        stop("vector or data.frame must be given")
    if (is.numeric(value)) {
        if (length(value) != nrow(object@design))
            stop("differing lengths")
        temp = data.frame(value)
        names(temp) = deparse(substitute(value))[1]
        value = temp
    }
    if (is.data.frame(value)) {
        if (nrow(value) != nrow(object@design))
            stop("differing number of rows")
    }
    object@response = value
    return(object)
})
setMethod(".nfp", "pbDesign", function(object) {
    x = factors(object)
    DB = FALSE
    if (is.list(x) && length(x[[1]]) > 0) {
        numAttr = length(attributes(x[[1]])) - 1
        .numFac = length(x)
        len = 0
        for (i in names(x)) if (length(x[[i]]@values) > len)
            len = length(x[[i]]@values)
        numAttr = numAttr + len
        numrows = numAttr - 1
        frameOut = data.frame(matrix(NA, ncol = .numFac, nrow = numrows))
        names(frameOut) = names(x)
        rownames(frameOut) = c(paste("value", 1:len), "name", "unit", "type")
        for (i in names(x)) {
            vin = 1:length(x[[i]]@values)
            frameOut[vin, i] = x[[i]]@values
            frameOut[numrows - 2, i] = x[[i]]@name
            frameOut[numrows - 1, i] = x[[i]]@unit
            frameOut[numrows, i] = x[[i]]@type
        }
        print(frameOut)
    }
})
setMethod("show", signature(object = "pbDesign"), function(object) {
    print(as.data.frame(object))
})
setMethod("summary", signature(object = "pbDesign"), function(object) {
    cat(paste("Plackett-Burman", toupper(object@designType), "Design"))
    cat("\n")
    cat("Information about the factors:\n\n")
    .nfp(object)
    cat("\n")
    cat("-----------\n")
    cat("\n")
    print(as.data.frame(object))
    cat("\n")
    cat("-----------\n")
    cat("\n")
})


.pbDesign=function(n)
{
  k=n-1
  if(k==27)
  {
  X=matrix(c(+1, -1, +1, +1, +1, +1, -1, -1, -1,
             +1, +1, -1, +1, +1, +1, -1, -1, -1,
             -1, +1, +1, +1, +1, +1, -1, -1, -1,
             -1, -1, -1, +1, -1, +1, +1, +1, +1,
             -1, -1, -1, +1, +1, -1, +1, +1, +1,
             -1, -1, -1, -1, +1, +1, +1, +1, +1,
             +1, +1, +1, -1, -1, -1, +1, -1, +1,
             +1, +1, +1, -1, -1, -1, +1, +1, -1,
             +1, +1, +1, -1, -1, -1, -1, +1, +1),nrow=9,ncol=9,byrow=TRUE)
  Y=matrix(c(-1, +1, -1, -1, -1, +1, -1, -1, +1,
             -1, -1, +1, +1, -1, -1, +1, -1, -1,
             +1, -1, -1, -1, +1, -1, -1, +1, -1,
             -1, -1, +1, -1, +1, -1, -1, -1, -1,
             +1, -1, -1, -1, -1, +1, +1, -1, -1,
             -1, +1, -1, +1, -1, -1, -1, +1, -1,
             -1, -1, +1, -1, -1, +1, -1, +1, -1,
             +1, -1, -1, +1, -1, -1, -1, -1, +1,
             -1, +1, -1, -1, +1, -1, +1, -1, -1),nrow=9,ncol=9,byrow=TRUE)
  Z=matrix(c(+1, +1, -1, +1, -1, +1, +1, -1, +1,
             -1, +1, +1, +1, +1, -1, +1, +1, -1,
             +1, -1, +1, -1, +1, +1, -1, +1, +1,
             +1, -1, +1, +1, +1, -1, +1, -1, +1,
             +1, +1, -1, -1, +1, +1, +1, +1, -1,
             -1, +1, +1, +1, -1, +1, -1, +1, +1,
             +1, -1, +1, +1, -1, +1, +1, +1, -1,
             +1, +1, -1, +1, +1, -1, -1, +1, +1,
             -1, +1, +1, -1, +1, +1, +1, -1, +1),nrow=9,ncol=9,byrow=TRUE)
  design=data.frame(rbind(cbind(X,Y,Z),cbind(X,Y,Z),cbind(X,Y,Z),rep(-1,27)))
  }
  else
  {
   if(k<3)
    stop("k needs to be grater than three!")
   if(k==3)
    firstRow=c(+1, -1, +1)[1:k]
   if(k>3 && k<=7)
    firstRow=c(+1, +1, +1, -1, +1, -1, -1)[1:k]
   if(k>=8 && k<=11)
    firstRow=c(+1, +1, -1, +1, +1, +1, -1, -1, -1, +1, -1)[1:k]
   if(k>=12 && k<=15)
    firstRow=c(+1, +1, +1, +1, -1, +1, -1, +1, +1, -1, -1, +1, -1, -1, -1)[1:k]
   if(k>=16 && k<=19)
    firstRow=c(+1, +1, -1, -1, +1, +1, +1, +1, -1, +1, -1, +1, -1, -1, -1, -1, +1, +1, -1)[1:k]
   if(k>=20 && k<=23)
    firstRow=c(+1, +1, +1, +1, +1, -1, +1, -1, +1, +1, -1, -1, +1, -1, -1, +1, -1, +1, -1, -1, -1, -1)[1:k]
   if(k>=24 && k<=26)
    firstRow=c(-1, -1, -1, -1, +1, -1, +1, -1, +1, +1, +1, -1, +1, +1, -1, -1, -1, +1, +1, +1, +1, +1, -1, -1, +1, +1, -1, +1, -1, -1, +1)[1:k]
   if(k==27)
    print("insert exception here!")
   if(k>=29 && k<=31)
    firstRow=c(-1, -1, -1, -1, +1, -1, +1, -1, +1, +1, +1, -1, +1, +1, -1, -1, -1, +1, +1, +1, +1, +1, -1, -1, +1, +1, -1, +1, -1, -1, +1)[1:k]
   if(k>=32 && k<=35)
    firstRow=c(-1, +1, -1, +1, +1, +1, -1, -1, -1, +1, +1, +1, +1, +1, -1, +1, +1, +1, -1, -1, +1, -1, -1, -1, -1, +1, -1, +1, -1, +1, +1, -1, -1, +1, -1)[1:k]
   if(k>=36 && k<=39)
    firstRow=c(+1, +1, -1, -1, +1, +1, +1, +1, -1, +1, -1, +1, -1, -1, -1, -1, +1, +1, -1, -1, +1, +1, -1, -1, +1, +1, +1, +1, -1, +1, -1, +1, -1, -1, -1, -1, +1, +1, -1)[1:k]
   if(k>=40 && k<=43)
    firstRow=c(+1, +1, -1, -1, +1, -1, +1, -1, -1, +1, +1, +1, -1, +1, +1, +1, +1, +1, -1, -1, -1, +1, -1, +1, +1, +1, -1, -1, -1, -1, -1, +1, -1, -1, -1, +1, +1, -1, +1, -1, +1, +1, -1)[1:k]
   if(k>=44 && k<=47)
    firstRow=c(+1, +1, +1, +1, +1, -1, +1, +1, +1, +1, -1, -1, +1, -1, +1, -1, +1, +1, +1, -1, -1, +1, -1, -1, +1, +1, -1, +1, +1, -1, -1, -1, +1, -1, +1, -1, +1, +1, -1, -1, -1, -1, +1, -1, -1, -1, -1)[1:k]
   if(k>=48 && k<=59)
    firstRow=c(+1, +1, -1, +1, +1, +1, -1, +1, -1, +1, -1, -1, +1, -1, -1, +1, +1, +1, -1, +1, +1, +1, +1, -1, -1, +1, +1, +1, +1, +1, -1, -1, -1, -1, -1, +1, +1, -1, -1, -1, -1, +1, -1, -1, -1, +1, +1, -1, +1, +1, -1, +1, -1, +1, -1, -1, -1, +1, -1)[1:k]

   design=matrix(NA,nrow=k+1,ncol=k)
   design[1,]=firstRow
   nextRow=firstRow

   for(i in 2:k)
   {
    nextRow=nextRow[c(k,1:k-1)]
    design[i,]=nextRow
   }
   lastRow=rep(-1,k)
   design[k+1,]=lastRow
   design=data.frame(design)
  }
 return(design)
}

pbDesign = function(n, k , randomize = TRUE, replicates = 1) {
 if(missing(n)&&missing(k))
  stop("Either n or k must be set!")
 if(missing(n)==FALSE && missing(k)==FALSE && k!=n-1 )
  stop("Wrong combination of n and k")
 if(missing(n))
  n=k+1
 if(missing(k))
  k=n-1  
 DB = FALSE
 odo = NA
 if (DB)
  print(n)
 design = .pbDesign(n)
 repVec = rep(1, nrow(design))
 if (replicates > 1) {
  X = .pbDesign(n)
   for (i in 1:(replicates - 1)) {
    design = rbind(design, X)
    repVec = c(repVec, rep(i + 1, times = nrow(X)))
   }
  }
 Replicate = data.frame(Replicate = as.numeric(repVec))
 if (DB)
  print(Replicate)
  odo = new("pbDesign")
  odo@design = design
#  names(odo@design) = .NAMES[1:ncol(design)]
  odo@replic = Replicate
  StandOrder = 1:nrow(odo@design)
  RunOrder = StandOrder
  if (randomize) {
   RunOrder = sample(1:nrow(odo@design), nrow(odo@design), replace = FALSE, prob = NULL)
  }
  odo@design = odo@design[order(RunOrder), ]
  odo@replic = data.frame(Replicate = odo@replic[order(RunOrder), 1])
  row.names(odo@design) = odo@design$RunOrder
  odo@runOrder = data.frame(RunOrder = data.frame(RunOrder = RunOrder)[order(RunOrder), ])
  odo@standardOrder = data.frame(StandOrder = data.frame(StandOrder = StandOrder)[order(RunOrder), ])
  odo@response = data.frame(y = rep(NA, nrow(odo@design)))
  tfList = vector("list", ncol(design))
  for (i in seq(along = tfList)) tfList[i] = new("pbFactor")
   names(tfList) = names(odo@design)
  factors(odo) = tfList
  valList = list(length = length(names(odo)))
  for (i in names(names(odo))) valList[[i]] = sort(unique(odo@design[, i]))
   values(odo) = valList
 return(odo)
}

