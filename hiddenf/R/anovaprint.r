anovaprint <-
function (x, digits = max(getOption("digits") - 2L, 3L), signif.stars = getOption("show.signif.stars"), 
    ...) 
{
    if (!is.null(heading <- attr(x, "heading"))) 
        cat(heading, sep = "\n")
    nc <- dim(x)[2L]
    if (is.null(cn <- colnames(x))) 
        stop("'anova' object must have colnames")
    has.P <- grepl("^(P|Pr)\\(", cn[nc])
    zap.i <- 1L:(if (has.P) 
        nc - 1
    else nc)
    i <- which(substr(cn, 2, 7) == " value")
    i <- c(i, which(!is.na(match(cn, c("F", "Cp", "Chisq")))))
    if (length(i)) 
        zap.i <- zap.i[!(zap.i %in% i)]
    tst.i <- i
    if (length(i <- grep("Df$", cn))) 
        zap.i <- zap.i[!(zap.i %in% i)]
    printCoefmat(x, digits = digits, signif.stars = signif.stars, 
        has.Pvalue = has.P, P.values = has.P, cs.ind = NULL, 
        zap.ind = zap.i, tst.ind = tst.i, na.print = "", ...)
    invisible(x)
}

