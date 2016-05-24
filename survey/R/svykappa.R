
svykappa<-function(formula, design,...) UseMethod("svykappa",design)

svykappa.default<-function(formula, design,...) {
    if (ncol(attr(terms(formula), "factors")) != 2) 
        stop("kappa is only computed for two variables")
    rows <- formula[[2]][[2]]
    cols <- formula[[2]][[3]]
    df <- model.frame(design)
    nrow <- length(unique(df[[as.character(rows)]]))
    ncol <- length(unique(df[[as.character(cols)]]))
    rnames<-paste(".",letters,"_",sep="")
    cnames<-paste(".",LETTERS,"_",sep="")
    if (nrow != ncol) 
        stop("number of categories is different")
    probs <- eval(bquote(svymean(~.(rows) + .(cols) + interaction(.(rows), 
        .(cols)), design)))
    nms <- c(rnames[1:nrow], cnames[1:ncol], outer(1:nrow, 
        1:ncol, function(i, j) paste(rnames[i], cnames[j], 
            sep = ".")))
    names(probs) <- nms
    v <- attr(probs, "var")
    dimnames(v) <- list(nms, nms)
    attr(probs, "var") <- v
    obs <- parse(text = paste(nms[nrow + ncol + 1+ (0:(nrow-1))*(ncol+1)], 
        collapse = "+"))[[1]]
    expect <- parse(text = paste(nms[1:nrow], nms[nrow + 1:ncol], 
        sep = "*", collapse = "+"))[[1]]
    svycontrast(probs, list(kappa = bquote((.(obs) - .(expect))/(1 - 
        .(expect)))))
}

