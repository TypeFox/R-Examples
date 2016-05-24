read.gct <-
function (file.gct, file.cls) 
{
    gct <- read.table(file.gct, skip = 2, sep = "\t", header = TRUE, 
        na.strings = "na", as.is = TRUE, quote = "\"")
    rownames(gct) <- gct[, 1]
    gct <- as.matrix(gct[, -(1:2)])
    cls <- readLines(file.cls)
    cls <- gsub("\t{1,}", " ", cls)
    cls <- strsplit(cls, split = " ")[[3]]
    names(cls) <- colnames(gct)
    gExprs <- list(gExprs = gct, sampleinfo = cbind(sampleid = colnames(gct), 
        status = cls))
    return(gExprs)
}
checkExprs <-
function (gExprs) 
{
    gexp <- gExprs$gExprs
    gexp <- t(apply(gexp, 1, function(Row) {
        if (length(ind.pos <- which(Row > 0)) == 0) {
            print("checkExprs(): one row has no positive values")
            return(rep(NA, length(Row)))
        }
        if (length(ind.na <- which(is.na(Row))) != 0) 
            Row[ind.na] <- mean(Row[ind.pos])
        posMin <- min(Row[ind.pos])
        if (length(ind.neg <- which(Row <= 0)) > 0) 
            Row[ind.neg] <- posMin
        return(Row)
    }))
    colnames(gexp) <- colnames(gExprs$gExprs)
    row.keep <- which(!is.na(gexp[, 1]))
    if (nrow.na <- nrow(gexp) - length(row.keep)) {
        gexp <- gexp[row.keep, ]
        print(paste(nrow.na, "rows are without positive values and removed."))
    }
    Med <- median(gexp)
    if (Med > 10) 
        gexp <- log2(gexp)
    gExprs$gExprs <- gexp
    return(gExprs)
}
collapseProbes <-
function (Mat, ann.data = NULL) 
{
    if (is.null(ann.data)) {
        ind.fil <- which(!is.na(rownames(Mat)))
        Mat <- Mat[ind.fil, ]
        row.factor <- as.factor(rownames(Mat))
    }
    else {
        ann.data <- ann.data[which(!is.na(ann.data))]
        ind.fil <- which(!is.na(match(rownames(Mat), names(ann.data))))
        Mat <- as.matrix(Mat[ind.fil, ])
        probe.ann <- ann.data[rownames(Mat)]
        row.factor <- as.factor(probe.ann)
    }
    exprs.collapsed <- tapply(1:nrow(Mat), row.factor, function(x) {
        if (length(x) > 1) 
            colMeans(Mat[x, ])
        else Mat[x, ]
    })
    rowNames <- names(exprs.collapsed)
    exprs.collapsed <- matrix(unlist(exprs.collapsed), byrow = T, 
        ncol = ncol(Mat))
    rownames(exprs.collapsed) <- rowNames
    colnames(exprs.collapsed) <- colnames(Mat)
    return(exprs.collapsed)
}
