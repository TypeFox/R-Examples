"boxes" <-
function (..., top = FALSE, shrink = 1, textcolor = NULL, yadj=NULL) 
{
    boxcall <- match.call()
    boxcall$top <- boxcall$shrink <- boxcall$textcolor <- boxcall$yadj <- NULL
    boxcall[[1]] <- as.name("boxplot")
    box <- eval(boxcall, parent.frame())
    mids <- 1:length(box$n)
    if (top) {
        where <- par("usr")[4]
        if (is.null(yadj)) {
            yadj <- 2
        }
        adj <- c(0.5, yadj)
    }
    else {
        where <- par("usr")[3]
        if (is.null(yadj)) {
            yadj <- -0.2
        }
        adj <- c(0.5, yadj)
    }
    tcex <- par("cex")
    par(cex = shrink * tcex)
    text(x = mids, y = where, labels = box$n, adj = adj, col = textcolor)
    par(cex = tcex)
    invisible(box)
}