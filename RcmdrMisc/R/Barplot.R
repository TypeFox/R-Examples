Barplot <- function(x, by, scale=c("frequency", "percent"), 
                    style=c("divided", "parallel"),
                    col=rainbow_hcl(length(levels(by))),
                    xlab=deparse(substitute(x)), 
                    legend.title=deparse(substitute(by)), ylab=scale,
                    legend.pos="topright", ...){
    if (!is.factor(x)) stop("x must be a factor")
    if (!missing(by) && !is.factor(by)) stop("by must be a factor")
    scale <- match.arg(scale)
    style <- match.arg(style)
    if (missing(by)){
        y <- table(x)
        if (scale == "percent") y <- 100*y/sum(y)
        barplot(y, xlab=xlab, ylab=ylab, ...)
    }
    else{
        col <- col[1:length(levels(by))]
        y <- table(by, x)
        if (scale == "percent") y <- 100*y/sum(y)
        barplot(y, xlab=xlab, ylab=ylab, legend.text=levels(by),
                col=col, 
                args.legend=list(x=legend.pos, title=legend.title, inset=0.05),
                beside = style == "parallel", ...)
    }
    return(invisible(NULL))
}
