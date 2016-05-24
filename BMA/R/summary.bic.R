# old (CF)
"print.bic.glm" <- 
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Posterior probabilities(%): \n")
    out <- x$probne0
    names(out) <- x$namesx
    print(out, ...)
    cat("\n Coefficient posterior expected values: \n")
    out <- x$postmean
    outnames<- c(NA, x$output.names)
    names(outnames)[1]<- "Intercept"
    nms <- NULL
    for (i in 1:length(outnames)) {
        if (is.na(outnames[i][1])) 
            nms <- c(nms, names(outnames[i]))
        else nms <- c(nms, paste(names(outnames[i]), unlist(outnames[i])[-1], 
            sep = "."))
    }
    names(out) <- nms
    
    fout<- format(out, digits=digits)
    fout[is.na(out)]<- ""
    print.default(fout, print.gap = 2, 
        quote = FALSE, ...)
    invisible(x)
}

print.bic.glm <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Posterior probabilities(%): \n")
    out <- x$probne0
#   names(out) <- if (x$factor.type) x$namesx else colnames(x$mle)[-1]
    print(out, ...)
    cat("\n Coefficient posterior expected values: \n")
    out <- x$postmean
#    names(out) <- c("(Intercept)", colnames(x$mle)[-1])
    fout <- format(out, digits = digits)
    fout[is.na(out)] <- ""
    print.default(fout, print.gap = 2, quote = FALSE, ...)
    invisible(x)
}

"print.bicreg" <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Posterior probabilities(%): \n")
    out <- x$probne0
    names(out) <- x$namesx
    print(out, ...)
    cat("\n Coefficient posterior expected values: \n")
    out <- x$postmean
    names(out) <- c("(Intercept)", x$namesx)
    print.default(format(out, digits = digits), print.gap = 2, 
        quote = FALSE, ...)
    invisible(x)
}

"print.bic.surv" <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("\n Posterior probabilities(%): \n")
    out <- x$probne0
    names(out) <- x$namesx
    print(out, ...)
    cat("\n Coefficient posterior expected values: \n")
    out <- x$postmean
    nms <- NULL
    for (i in 1:length(x$output.names)) {
        if (is.na(x$output.names[i][1])) 
            nms <- c(nms, names(x$output.names[i]))
        else nms <- c(nms, paste(names(x$output.names[i]), unlist(x$output.names[i])[-1], 
            sep = "."))
    }
    names(out) <- nms
    print.default(format(out, digits = digits), print.gap = 2, 
        quote = FALSE, ...)
    invisible(x)
}




"summary.bic.glm" <-
function (object, n.models = 5, digits = max(3, getOption("digits") - 
    3), conditional = FALSE, display.dropped = FALSE, ...) 
{
    x<- object
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (display.dropped & x$reduced) {
        cat("\nThe following variables were dropped prior to averaging:\n")
        cat(x$dropped)
        cat("\n")
    }
    n.models <- min(n.models, x$n.models)
    sel <- 1:n.models
    cat("\n ", length(x$postprob), " models were selected")
    cat("\n Best ", n.models, " models (cumulative posterior probability = ", 
        round(sum(x$postprob[sel]), digits), "): \n\n")

    x$namesx<- c("Intercept", x$namesx)
    nms <- length(x$namesx)
    ncx <- length(unlist(x$assign))
    nvar <- rep(0, times = n.models)
    for (i in 1:(nms-1)) nvar <- nvar + as.numeric(as.vector(rbind(rep(1, 
        length(x$assign[[i+1]]))) %*% (t(x$mle[sel, x$assign[[i+1]], drop = FALSE] != 0)) > 0))



    modelposts <- format(round(x$postprob[sel], 3), digits = 3)
    coeffs <- t(x$mle[sel, , drop = FALSE])
    cfbic <- rbind(x$bic[sel], coeffs)
    cfbicf <- format(cfbic, digits = digits)
    coeffsf <- cfbicf[-1, , drop = FALSE]
    bic <- cfbicf[1, , drop = FALSE]
    postmeans <- format(x$postmean, digits = digits)
    postsds <- format(x$postsd, digits = digits)

    postmeans[is.na(x$postmean)]<- ""
    postsds[is.na(x$postsd)]<- ""

    if (conditional) {
        cpostmeans <- format(x$condpostmean, digits = digits)
        cpostsds <- format(x$condpostsd, digits = digits)
        cpostmeans[is.na(x$condpostmean)]<- ""
        cpostsds[is.na(x$condpostsd)]<- ""

    }
    varposts <- format(round(x$probne0, 1), digits = 3)
    strlength <- nchar(coeffsf[1, 1])
    decpos <- nchar(unlist(strsplit(coeffsf[2, 1], "\\."))[1])
    offset <- paste(rep(" ", times = decpos - 1), sep = "", collapse = "")
    offset2 <- paste(rep(" ", times = decpos + 1), sep = "", 
        collapse = "")
    modelposts <- paste(offset, modelposts, sep = "")
    nvar <- paste(offset2, nvar, sep = "")
    dotoffset <- round(max(nchar(coeffsf))/2)
    zerocoefstring <- paste(paste(rep(" ", times = dotoffset), 
        collapse = "", sep = ""), ".", sep = "")
    coeffsf[coeffs == 0] <- zerocoefstring
    coeffsf[is.na(coeffs)]<- ""
    avp <- NULL

    outnames<- c(NA, x$output.names)
    names(outnames)[1]<- "Intercept"
    varposts<- c("100",varposts)

    for (i in 1:nms) {
        avp <- rbind(avp, varposts[i])
        if (!is.na(outnames[[i]][1])) 
            avp <- rbind(avp, cbind(rep("", times = length(x$assign[[i]]))))
    }
    top <- cbind(postmeans, postsds)
    if (conditional) 
        top <- cbind(top, cpostmeans, cpostsds)
    top <- cbind(top, coeffsf)
    atop <- NULL
    for (i in 1:nms) {
        if (!is.na(outnames[[i]][1])) 
            atop <- rbind(atop, rbind(rep("", times = ncol(top))))
        atop <- rbind(atop, top[x$assign[[i ]], ])
    }
    top <- cbind(avp, atop)
    linesep <- rep("", times = ncol(top))
    offset <- c("", "", "")
    if (conditional) 
        offset <- c(offset, c("", ""))
    bottom <- rbind(c(offset, nvar), c(offset, bic), c(offset, 
        modelposts))
    out <- rbind(top, linesep, bottom)
    vnames <- NULL
    for (i in 1:nms) {
        vnames <- c(vnames, names(outnames[i]))
        blnk <- paste(rep(" ", times = nchar(names(outnames[i]))), 
            collapse = "")
        if (!is.na(outnames[i][1])) 
            vnames <- c(vnames, paste(blnk, unlist(outnames[i])[-1], 
                sep = "."))
    }
    row.names(out) <- c(vnames, "", "nVar", "BIC", "post prob")
    colnms <- c("p!=0", " EV", "SD")
    if (conditional) 
        colnms <- c(colnms, "cond EV", "cond SD")
    colnms <- c(colnms, paste("model ", 1:n.models, sep = ""))
    dimnames(out)[[2]] <- colnms
    print.default(out, print.gap = 2, quote = FALSE, ...)
}




"summary.bicreg" <-
function (object, n.models = 5, digits = max(3, getOption("digits") - 
    3), conditional = FALSE, display.dropped = FALSE, ...) 
{
    x<- object
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (display.dropped & x$reduced) {
        cat("\nThe following variables were dropped prior to averaging:\n")
        cat(x$dropped)
        cat("\n")
    }
    n.models <- min(n.models, x$n.models)
    sel <- 1:n.models
    cat("\n ", length(x$postprob), " models were selected")
    cat("\n Best ", n.models, " models (cumulative posterior probability = ", 
        round(sum(x$postprob[sel]), digits), "): \n\n")
    nms <- length(x$namesx) + 1
    r2 <- format(round(x$r2[sel]/100, 3), digits = 3)
    nvar <- rbind(rep(1, length(x$namesx) + 1)) %*% t(x$ols[sel, 
        ] != 0) - 1
    modelposts <- format(round(x$postprob[sel], 3), digits = 3)
    coeffs <- t(x$ols[sel, ])
    cfbic <- rbind(x$bic[sel], coeffs)
    cfbicf <- format(cfbic, digits = digits)
    coeffsf <- cfbicf[-1, ]
    bic <- cfbicf[1, ]
    dotoffset <- round(max(nchar(coeffsf))/2)
    zerocoefstring <- paste(paste(rep(" ", times = dotoffset), 
        collapse = "", sep = ""), ".", sep = "")
    coeffsf[coeffs == 0] <- zerocoefstring
    postmeans <- format(x$postmean, digits = digits)
    postsds <- format(x$postsd, digits = digits)
    if (conditional) {
        cpostmeans <- format(x$condpostmean, digits = digits)
        cpostsds <- format(x$condpostsd, digits = digits)
    }
    varposts <- format(round(c(100, x$probne0), 1), digits = 3)
    strlength <- nchar(coeffsf[1, 1])
    decpos <- nchar(unlist(strsplit(coeffsf[1, 1], "\\."))[1])
    offset <- paste(rep(" ", times = decpos - 1), sep = "", collapse = "")
    offset2 <- paste(rep(" ", times = decpos + 1), sep = "", 
        collapse = "")
    r2 <- paste(offset, r2, sep = "")
    modelposts <- paste(offset, modelposts, sep = "")
    nvar <- paste(offset2, nvar, sep = "")
    top <- cbind(varposts, postmeans, postsds)
    if (conditional) 
        top <- cbind(top, cpostmeans, cpostsds)
    top <- cbind(top, coeffsf)
    linesep <- rep("", times = ncol(top))
    offset <- c("", "", "")
    if (conditional) 
        offset <- c(offset, c("", ""))
    bottom <- rbind(c(offset, nvar), c(offset, r2), c(offset, 
        bic), c(offset, modelposts))
    out <- rbind(top, linesep, bottom)
    row.names(out) <- c("Intercept", x$namesx, "", "nVar", 
        "r2", "BIC", "post prob")
    colnms <- c("p!=0", " EV", "SD")
    if (conditional) 
        colnms <- c(colnms, "cond EV", "cond SD")
    colnms <- c(colnms, paste("model ", 1:n.models, sep = ""))
    dimnames(out)[[2]] <- colnms
    print.default(out, print.gap = 2, quote = FALSE, ...)
}
"summary.bic.surv" <-
function (object, n.models = 5, digits = max(3, getOption("digits") - 
    3), conditional = FALSE, display.dropped = FALSE, ...) 
{
    x<- object
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (display.dropped & x$reduced) {
        cat("\nThe following variables were dropped prior to averaging:\n")
        cat(x$dropped)
        cat("\n")
    }
    n.models <- min(n.models, x$n.models)
    sel <- 1:n.models
    cat("\n ", length(x$postprob), " models were selected")
    cat("\n Best ", n.models, " models (cumulative posterior probability = ", 
        round(sum(x$postprob[sel]), digits), "): \n\n")
    nms <- length(x$namesx)
    ncx <- length(unlist(x$assign)[-1])
    nvar <- rep(0, times = n.models)
    for (i in 1:nms) nvar <- nvar + as.numeric(as.vector(rbind(rep(1, 
        length(x$assign[[i + 1]]))) %*% (t(x$mle[sel, x$assign[[i + 
        1]], drop = FALSE] != 0)) > 0))
    modelposts <- format(round(x$postprob[sel], 3), digits = 3)
    coeffs <- t(x$mle[sel, , drop = FALSE])
    cfbic <- rbind(x$bic[sel], coeffs)
    cfbicf <- format(cfbic, digits = digits)
    coeffsf <- cfbicf[-1, , drop = FALSE]
    bic <- cfbicf[1, , drop = FALSE]
    postmeans <- format(x$postmean, digits = digits)
    postsds <- format(x$postsd, digits = digits)
    if (conditional) {
        cpostmeans <- format(x$condpostmean, digits = digits)
        cpostsds <- format(x$condpostsd, digits = digits)
    }
    varposts <- format(round(x$probne0, 1), digits = 3)
    strlength <- nchar(coeffsf[1, 1])
    decpos <- nchar(unlist(strsplit(coeffsf[1, 1], "\\."))[1])
    offset <- paste(rep(" ", times = decpos - 1), sep = "", collapse = "")
    offset2 <- paste(rep(" ", times = decpos + 1), sep = "", 
        collapse = "")
    modelposts <- paste(offset, modelposts, sep = "")
    nvar <- paste(offset2, nvar, sep = "")
    dotoffset <- round(max(nchar(coeffsf))/2)
    zerocoefstring <- paste(paste(rep(" ", times = dotoffset), 
        collapse = "", sep = ""), ".", sep = "")
    coeffsf[coeffs == 0] <- zerocoefstring
    avp <- NULL
    for (i in 1:nms) {
        avp <- rbind(avp, varposts[i])
        if (!is.na(x$output.names[[i]][1])) 
            avp <- rbind(avp, cbind(rep("", times = length(x$assign[[i + 
                1]]))))
    }
    top <- cbind(postmeans, postsds)
    if (conditional) 
        top <- cbind(top, cpostmeans, cpostsds)
    top <- cbind(top, coeffsf)
    atop <- NULL
    for (i in 1:nms) {
        if (!is.na(x$output.names[[i]][1])) 
            atop <- rbind(atop, rbind(rep("", times = ncol(top))))
        atop <- rbind(atop, top[x$assign[[i + 1]], ])
    }
    top <- cbind(avp, atop)
    linesep <- rep("", times = ncol(top))
    offset <- c("", "", "")
    if (conditional) 
        offset <- c(offset, c("", ""))
    bottom <- rbind(c(offset, nvar), c(offset, bic), c(offset, 
        modelposts))
    out <- rbind(top, linesep, bottom)
    vnames <- NULL
    for (i in 1:nms) {
        vnames <- c(vnames, names(x$output.names[i]))
        blnk <- paste(rep(" ", times = nchar(names(x$output.names[i]))), 
            collapse = "")
        if (!is.na(x$output.names[i][1])) 
            vnames <- c(vnames, paste(blnk, unlist(x$output.names[i])[-1], 
                sep = "."))
    }
    row.names(out) <- c(vnames, "", "nVar", "BIC", "post prob")
    colnms <- c("p!=0", " EV", "SD")
    if (conditional) 
        colnms <- c(colnms, "cond EV", "cond SD")
    colnms <- c(colnms, paste("model ", 1:n.models, sep = ""))
    dimnames(out)[[2]] <- colnms
    print.default(out, print.gap = 2, quote = FALSE, ...)
}

# duplicate (CF ?)
`summary.bicreg` <-
function (object, n.models = 5, digits = max(3, getOption("digits") - 
    3), conditional = FALSE, display.dropped = FALSE, ...) 
{
    x <- object
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (display.dropped & x$reduced) {
        cat("\nThe following variables were dropped prior to averaging:\n")
        cat(x$dropped)
        cat("\n")
    }
    n.models <- min(n.models, x$n.models)
    sel <- 1:n.models
    cat("\n ", length(x$postprob), " models were selected")
    cat("\n Best ", n.models, " models (cumulative posterior probability = ", 
        round(sum(x$postprob[sel]), digits), "): \n\n")
    nms <- length(x$namesx) + 1
    r2 <- format(round(x$r2[sel]/100, 3), digits = 3)
    nvar <- rbind(rep(1, length(x$namesx) + 1)) %*% t(x$ols[sel, , drop = FALSE
        ] != 0) - 1
    modelposts <- format(round(x$postprob[sel], 3), digits = 3)
    coeffs <- t(x$ols[sel, ,drop=FALSE])
    cfbic <- rbind(x$bic[sel], coeffs)
    cfbicf <- format(cfbic, digits = digits)
    coeffsf <- cfbicf[-1, ,drop=FALSE]
    bic <- cfbicf[1, ]
    dotoffset <- round(max(nchar(coeffsf))/2)
    zerocoefstring <- paste(paste(rep(" ", times = dotoffset), 
        collapse = "", sep = ""), ".", sep = "")
    coeffsf[coeffs == 0] <- zerocoefstring
    postmeans <- format(x$postmean, digits = digits)
    postsds <- format(x$postsd, digits = digits)
    if (conditional) {
        cpostmeans <- format(x$condpostmean, digits = digits)
        cpostsds <- format(x$condpostsd, digits = digits)
    }
    varposts <- format(round(c(100, x$probne0), 1), digits = 3)
    strlength <- nchar(coeffsf[1, 1])
    decpos <- nchar(unlist(strsplit(coeffsf[1, 1], "\\."))[1])
    offset <- paste(rep(" ", times = decpos - 1), sep = "", collapse = "")
    offset2 <- paste(rep(" ", times = decpos + 1), sep = "", 
        collapse = "")
    r2 <- paste(offset, r2, sep = "")
    modelposts <- paste(offset, modelposts, sep = "")
    nvar <- paste(offset2, nvar, sep = "")
    top <- cbind(varposts, postmeans, postsds)
    if (conditional) 
        top <- cbind(top, cpostmeans, cpostsds)
    top <- cbind(top, coeffsf)
    linesep <- rep("", times = ncol(top))
    offset <- c("", "", "")
    if (conditional) 
        offset <- c(offset, c("", ""))
    bottom <- rbind(c(offset, nvar), c(offset, r2), c(offset, 
        bic), c(offset, modelposts))
    out <- rbind(top, linesep, bottom)
    row.names(out) <- c("Intercept", x$namesx, "", "nVar", "r2", 
        "BIC", "post prob")
    colnms <- c("p!=0", " EV", "SD")
    if (conditional) 
        colnms <- c(colnms, "cond EV", "cond SD")
    colnms <- c(colnms, paste("model ", 1:n.models, sep = ""))
    dimnames(out)[[2]] <- colnms
    print.default(out, print.gap = 2, quote = FALSE, ...)
}

