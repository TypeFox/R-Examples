HGLMFactorList <-
function (formula, fr, rmInt, drop) 
{
    mf <- fr$mf
    bars <- expandSlash(findbars(formula[[3]]))
    for (i in 1:length(bars)) {
        checkcorr <- findplus(bars[[i]])
        if (checkcorr == 1) 
            stop("Correlated random effects are not currently allowed in the HGLM routines")
        if (checkcorr == -1) 
            stop("You do not need to specify '-1' for no intercept it is done be default")
    }
    if (!length(bars)) 
        stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars, function(x) {
        ff <- eval(substitute(as.factor(fac)[, drop = TRUE], 
            list(fac = x[[3]])), mf)
        im <- as(ff, "sparseMatrix")
        if (!isTRUE(validObject(im, test = TRUE))) 
            stop("invalid conditioning factor in random effect: ", 
                format(x[[3]]))
        if (is.name(x[[2]])) {
            tempexp <- paste("~", as.character(x[[2]]), "-1")
            tempexp <- as.formula(tempexp)[[2]]
        }
        else tempexp <- x[[2]]
        mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), 
            mf)
        if (rmInt) {
            if (is.na(icol <- match("(Intercept)", colnames(mm)))) 
                break
            if (ncol(mm) < 2) 
                stop("lhs of a random-effects term cannot be an intercept only")
            mm <- mm[, -icol, drop = FALSE]
        }
        ans <- list(f = ff, A = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) im)), Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) {
                im@x <- mm[, j]
                im
            })), ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
            colnames(mm))))
        if (drop) {
            ans$A@x <- rep(0, length(ans$A@x))
            ans$Zt <- drop0(ans$Zt)
        }
        ans
    })
    Design <- list(0)
    Subject <- list(0)
    for (i in 1:length(fl)) {
        Subject[[i]] <- as.factor(fl[[i]]$f)
        tempmat <- fl[[i]]$Zt
        tempmat <- as.matrix(t(tempmat))
        Design[[i]] <- tempmat
    }
    list(Design = Design, Subject = Subject, namesRE = names(bars))
}

