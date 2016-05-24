`image.mefa` <-
function(x, segm=NULL, trafo=c("none", "log", "bins", "prab"), 
probs = seq(0, 1, 0.05), ordering=TRUE, reverse=TRUE, names = FALSE,
show=TRUE, ylab, xlab, ...)
{
if (!is.mefa(x))
    stop("object is not of class 'mefa'")
if (length(trafo) > 1) trafo <- trafo[1]
trafo <- match.arg(trafo, c("none", "log", "bins", "prab"))
m <- if (!is.null(segm))
    x$segm[[segm]] else x$xtab
# ordering based on totals
if (ordering)
    m <- m[order(rowSums(x$xtab)), order(colSums(x$xtab),decreasing=TRUE)]
mm <- t(m)
# transformations
if (trafo == "log")
    mm <- log10(mm)
if (trafo == "bins")
    mm <- matrix(qvector(array(mm), probs), nrow(mm), ncol(mm))
if (trafo == "prab")
    mm <- matrix(as.numeric(mm > 0), nrow(mm), ncol(mm))
if (reverse)
    mm <- max(mm) - mm

if (missing(ylab)) ylab <- "Samples"
if (missing(xlab)) xlab <- "Taxa"

if (show) {
    image(1:ncol(m), 1:nrow(m), mm, ylab=ylab, xlab=xlab, axes = FALSE, ...)
    box()
    if (length(names) == 1)
        names <- rep(names, 2)
    if (names[1]) {
        axis(2, at = 1:ncol(mm), labels = colnames(mm),
            las = 2, ...)
    }
    if (names[2]) {
        axis(1, at = 1:nrow(mm), labels = rownames(mm),
            las = 2, ...)
    }
}
if (show)
    invisible(mm) else return(mm)
}
