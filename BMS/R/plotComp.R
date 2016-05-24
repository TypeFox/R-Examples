plotComp <-
function (..., varNr = NULL, comp = "PIP", exact = FALSE, include.legend = TRUE, 
    add.grid = TRUE, do.par = TRUE, cex.xaxis = 0.8) 
{
    col_default = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
        "#66A61E", "#E6AB02", "#A6761D", "#666666")
    bmaList = list(...)
    bmaix = sapply(bmaList, is.bma)
    bmaList = bmaList[bmaix]
    bmaNr = length(bmaList)
    if (!all(sapply(bmaList, is.bma))) {
        stop("Submit only bma objects to compare results (no other objects)")
    }
    dotargs = match.call(expand.dots = FALSE)$...
    dotargs = .adjustdots(dotargs, ylab = paste(comp), pch = 1:bmaNr, 
        col = col_default, type = "p", lty = 1:5, lwd = 1.5, 
        xlab = "", xaxt = "n")
    dotargs = dotargs[!c(bmaix, logical(length(dotargs) - length(bmaix)))]
    xMat = lapply(bmaList, function(x) rownames(estimates.bma(x, 
        exact = exact)))
    xNames = xMat[[1]]
    ind = as.numeric(unlist(lapply(xMat, function(x) length(x))))
    if (length(unique(ind) > 1)) {
        smallestSet = which.min(ind)
        indMat = array(0:0, dim = c(length(xMat[[smallestSet]]), 
            length(ind)))
        for (i in 1:length(ind)) {
            indMat[, i] = as.numeric(xMat[[smallestSet]] %in% 
                xMat[[i]])
        }
        xNamesInd = which(rowSums(indMat) == bmaNr)
        xNames = xMat[[smallestSet]][xNamesInd]
    }
    compNames = c(colnames(estimates.bma(bmaList[[1]])), "Std Mean", 
        "Std Coef")
    if (is.null(xNames)) {
        stop("the bma objects have to have (the same) rownames attached to them")
    }
    if (!(comp %in% compNames)) {
        stop("Please specify comp as one of PIP, Post Mean, Post SD, Std Mean, or Std Coef")
    }
    if (comp == "Std Mean") {
        compMatrix = sapply(bmaList, function(x) estimates.bma(x, 
            std.coefs = TRUE, exact = exact)[xNames, "Post Mean"])
        comp = "Standardized Coefficients"
    }
    else if (comp == "Std SD") {
        compMatrix = sapply(bmaList, function(x) estimates.bma(x, 
            std.coefs = TRUE, exact = exact)[xNames, "Post SD"])
        comp = "Standardized SD"
    }
    else {
        compMatrix = sapply(bmaList, function(x) estimates.bma(x, 
            exact = exact)[xNames, comp])
    }
    bmaNames = names(list(...))[bmaix]
    colnames(compMatrix) = paste("Model", 1:bmaNr)
    if (!is.null(bmaNames) && (length(bmaNames) == ncol(compMatrix))) {
        for (bix in 1:bmaNr) {
            colnames(compMatrix)[[bix]] <- ifelse(bmaNames[[bix]] == 
                "", paste("Model", bix), bmaNames[[bix]])
        }
    }
    if (!is.null(varNr)) {
        compMatrix = compMatrix[varNr, , drop = FALSE]
    }
    if (as.logical(do.par)) {
        oldmar = par()$mar
        spaceforxaxis = strwidth(rownames(compMatrix)[which.max(nchar(rownames(compMatrix)))], 
            units = "inches", cex = cex.xaxis) * (par("mar")/par("mai"))[[2]]
        tempmar = oldmar
        tempmar[1] = min(max(oldmar[1], spaceforxaxis + oldmar[1]/3), 
            0.5 * par("fin")[[2]] * (par("mar")/par("mai"))[[1]])
        par(mar = tempmar)
    }
    eval(as.call(c(list(as.name("matplot"), as.name("compMatrix")), 
        as.list(dotargs))))
    if (as.logical(include.legend)) {
        extractfromdotargs = function(...) {
            dal = list(...)
            return(list(col = dal$col, pch = dal$pch))
        }
        myargs = eval(as.call(c(list(as.name("extractfromdotargs")), 
            as.list(dotargs))))
        legend("topright", colnames(compMatrix), pch = myargs$pch, 
            col = myargs$col, bty = "n")
    }
    if (as.logical(add.grid)) 
        grid()
    axis(1, las = 2, at = 1:nrow(compMatrix), labels = rownames(compMatrix), 
        cex.axis = cex.xaxis)
    if (as.logical(do.par)) 
        par(mar = oldmar)
}
