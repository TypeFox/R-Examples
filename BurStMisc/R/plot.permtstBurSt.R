"plot.permtstBurSt" <-
function (x, col = c("black", "red"), width = 10, uniqlim = 10, 
    main = "", xlab = "Scores", ...) 
{
    orig.score <- x$original.score
    ulen <- length(unique(x$perm.scores))
    if (ulen > uniqlim) {
        hist(x$perm.scores, xlim = range(x$perm.scores, orig.score), 
            main = main, xlab = xlab, ...)
        box()
        abline(v = orig.score, col = col[2])
    }
    else {
        ptab <- table(x$perm.scores)
        vals <- as.numeric(names(ptab))
        if (x$alternative == "greater") {
            extreme <- vals >= orig.score
        }
        else {
            extreme <- vals <= orig.score
        }
        if (all(vals > orig.score) || all(vals < orig.score)) {
            xrng <- range(vals, orig.score)
        }
        else {
            xrng <- range(vals)
        }
        plot(vals, ptab, type = "n", xlab = xlab, ylab = "Count", 
            xlim = xrng, ...)
        points(vals[!extreme], ptab[!extreme], type = "h", col = col[1], 
            lwd = width, ...)
        points(vals[extreme], ptab[extreme], type = "h", col = col[2], 
            lwd = width, ...)
        if (nchar(main)) 
            title(main = main)
    }
}
