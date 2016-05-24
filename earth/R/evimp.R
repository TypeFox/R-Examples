# evimp.R: estimate variable importances in an earth object

# Return a vector of column numbers for predictors that are used
# in the final model

get.used.preds <- function(object)   # object is an earth object
{
    which(apply(object$dirs[object$selected.terms,,drop=FALSE],2,any1))
}
# Print predictors in order of decreasing estimated importance.
# A one line summary.  Called by print.summary.earth.

print.one.line.evimp <- function(object) # object is an "earth" object
{
    if(is.null(object$prune.terms)) {
        if(is.null(object$ifold)) { # not a fold model?
            # must have been created by mars.to.earth
            cat("Importance: object has no prune.terms, call update() on the model to fix that\n")
        }
        return()
    }
    evimp <- row.names(evimp(object, trim=FALSE))
    if(length(evimp) == 0)
        cat("Importance: no predictors")
    else if(length(evimp) == 1)
        cat0("Importance: ", evimp[1])
    else {
        width <- max(getOption("width")-5, 20) # -5 for ", ..."
        s <- paste0("Importance: ", evimp[1])
        for(ipred in 2:length(evimp)) {
            temp <- paste0(s, ", ", evimp[ipred])
            if(nchar(temp) >= width) {
                s <- paste0(s, ", ...")
                break
            }
            s <- temp
        }
        cat(s)
    }
    cat("\n")
}
evimp <- function(object, trim=TRUE, sqrt.=TRUE) # see help page for description
{
    trim  <- check.boolean(trim)
    sqrt. <- check.boolean(sqrt.)

    # convert col numbers in predtab to col numbers in importances
    as.icriti <- function(icrit) c(3,4,6)[icrit]

    check.classname(object, substitute(object), "earth")
    stopifnot(!is.null(object$prune.terms))
    nsubsets <- length(object$selected.terms)
    dirs <- object$dirs
    pred.names <- gen.colnames(dirs, "x", "x", trace=0)

    # tagged.pred.names is a copy of pred.names but with unused
    # predictors renamed by adding a "-unused" suffix.
    # By unused, we mean unused in the final model.

    used.preds <- to.logical(get.used.preds(object), len=length(pred.names))
    tagged.pred.names <- pred.names
    tagged.pred.names[!used.preds] <-
            paste0(tagged.pred.names[!used.preds], "-unused")

    # deltas[isubset, icrit] is the change in criterion value
    # for isubset using criterion icrit

    stopifnot(nsubsets >= 1)
    deltas <- matrix(nrow=nsubsets-1, ncol=3)
    colnames(deltas) <- c("nsubsets", "gcv", "rss")
    deltas[,"nsubsets"] <- rep(1, times=nsubsets-1)
    deltas[,"gcv"]      <- -diff(object$gcv.per.subset[seq_len(nsubsets)])
    deltas[,"rss"]      <- -diff(object$rss.per.subset[seq_len(nsubsets)])

    # preds.in.each.term[iterm] is the indices of predictors in term iterm

    preds.in.each.term <- apply(object$dirs, 1, function(row) which(row != 0))

    # importances is the matrix we return

    importances <- matrix(0, nrow=length(pred.names), ncol=7)
    colnames(importances) <- c("col", "used", "nsubsets", "gcv", "gcv.match", "rss", "rss.match")
    rownames(importances) <- tagged.pred.names
    importances[, "col"] <- seq_len(nrow(importances))
    importances[used.preds, "used"] <- 1

    if(nsubsets > 1) for(isubset in 2:nsubsets) {
        terms.in.this.subset <- object$prune.terms[isubset,-1]  # -1 drops intercept
        preds.in.this.subset <-
            unique(unlist(preds.in.each.term[terms.in.this.subset]))

        for(icrit in 1:3) {
            icriti <- as.icriti(icrit)
            importances[preds.in.this.subset, icriti] <-
                importances[preds.in.this.subset, icriti] +
                deltas[isubset-1, icrit]
        }
    }
    # sort rows in "importances" by the nsubsets criteria
    # and with the "gcv" criterion as a secondary sort key

    order.nsubsets <- order(importances[,"nsubsets"], importances[,"gcv"], decreasing=TRUE)
    importances <- importances[order.nsubsets, , drop=FALSE]

    if(nrow(importances) > 1)
        for(icrit in 2:3) {
            # tag importances where gcv or rss ordering disagrees with nsubsets ordering

            icriti <- as.icriti(icrit)
            importances[, icriti+1] <- 1
            for(i in 2:nrow(importances))
                if(importances[i,icriti] > importances[i-1,icriti])
                    importances[i, icriti+1] <- 0

            # normalize importances

            max <- max(abs(importances[,icriti]))
            if(max != 0) {
                if(sqrt.) {
                    temp <- sqrt(abs(importances[,icriti]) / max)
                    signs <- ifelse(importances[,icriti] < 0, -1, 1)
                    importances[,icriti] <- 100 * signs * temp
                } else
                    importances[,icriti] <- 100 * importances[,icriti] / max
            }
        }

    if(trim) {
        # keep only rows for predictors that are used in at least one subset

        in.at.least.one.subset <- importances[,"nsubsets"] != 0
        importances <- importances[in.at.least.one.subset, , drop=FALSE]
    }
    class(importances) <- "evimp"   # allows use of plot.evimp
    attr(importances, "sqrt") <- sqrt.
    importances
}
print.evimp <- function(x = stop("no 'x' argument"), ...) # x is an "evimp" object
{
    stopifnot(NCOL(x) == 7)
    if(NROW(x) == 0) {
        printf("    nsubsets   gcv    rss\n")
        return()
    }
    # truncate rownames if necessary so each entry requires only one line on the screen
    rownames <- rownames(x)
    max.rowname <- max(nchar(rownames))
    width <- getOption("width")
    if(max.rowname > width-25)  { # width of stuff to right of rowname is slighty less than 25
        rownames <- substr(rownames, 1, max(20, width-25))
        max.rowname <- max(nchar(rownames))
    }
    printf("%*s nsubsets   gcv    rss\n", max.rowname, " ")
    for(i in seq_len(nrow(x)))
        printf("%-*s %8d %5.1f%s %5.1f%s\n",
            max.rowname, rownames[i], x[i, 3],
            x[i, 4], if(x[i, 7]) " " else ">",
            x[i, 6], if(x[i, 7]) "" else ">")
}
# TODO this would be better if rotated clockwise 90 degrees so could easily read var names

plot.evimp <- function(
    x = stop("no 'x' argument"),
    cex.var = 1,

    type.nsubsets = "l",
    col.nsubsets = "black",
    lty.nsubsets = 1,

    type.gcv = "l",
    col.gcv = 2,
    lty.gcv = 1,

    type.rss = "l",
    col.rss = "gray60",
    lty.rss = 1,

    cex.legend = 1,
    x.legend = nrow(x),
    y.legend = x[1,"nsubsets"],

    rh.col = 1,
    do.par = TRUE,
    ...)
{
    check.classname(x, substitute(x), "evimp")
    # make sure that all evimp columns are present (extra columns are ok)
    if(any(pmatch(c("col", "used", "nsubsets", "gcv"), colnames(x), nomatch=0) == 0))
        stop0("x is not an evimp matrix")
    if(nrow(x) == 0) { # intercept-only model
        max.subsets <- 0
        varlabs <- "intercept"
    } else {
        max.subsets <- x[1, "nsubsets"]
        varlabs <- paste(rownames(x), sprintf("%3d", x[,"col"]))
    }
    sqrt. <- if(attr(x, "sqrt", exact=TRUE)) TRUE else FALSE
    par <- par("mar", "cex")
    on.exit(par(par))
    cex.var <- par$cex * cex.var    # cex.var is relative to current cex
    do.par <- check.boolean(do.par)
    if(do.par) {
        # TODO what is the best way of doing the bottom.margin calculation?
        # The .5 is a hack to convert nchars to line heights, as required by mar
        mar <- par$mar
        mar[1] <- cex.var * .5 * max(nchar(varlabs) + 6)    # bottom margin
        mar[4] <- mar[4] + 3                                # right margin
        par(mar=mar) # big bottom and right margins
    }
    main <- dot("main", DEF="Variable importance", ...)
    if(max.subsets == 0) {
        plot(1, ylim=c(0, 1), type=type.nsubsets, # intercept-only model, dummy plot
             xlab="", xaxt="n", ylab="nsubsets",
             main=main, lty=lty.nsubsets, col=col.nsubsets)
    } else {
        plot(x[, "nsubsets"], ylim=c(0, max.subsets), type=type.nsubsets,
             xlab="", xaxt="n", ylab="nsubsets",
             main=main, lty=lty.nsubsets, col=col.nsubsets)
        lines(max.subsets * x[,"rss"] / 100, type=type.rss, lty=lty.rss, col=col.rss)
        # plot gcv second so it goes on top of rss (gcv arguably more important than rss)
        lines(max.subsets * x[,"gcv"] / 100, type=type.gcv, lty=lty.gcv, col=col.gcv)
    }
    zero.or.one.var <- nrow(x) <= 1
    if(is.specified(x.legend)) {
        if(sqrt.)
            legend <- c("nsubsets", "sqrt gcv", "sqrt rss")
        else
            legend <- c("nsubsets", "gcv", "rss")
        legend(x=if(zero.or.one.var) "topright" else x.legend,
               y = y.legend, xjust=1,
               legend=legend,
               col=c(col.nsubsets, col.gcv, col.rss),
               lty=c(lty.nsubsets, lty.gcv, lty.rss),
               bg="white", cex=cex.legend)
    }
    # right hand axis: normalized rss/gcv values, always 0...100
    # TODO how to get the x position in the call to text correct for all window sizes?
    axis(side=4,
         at=c(0,.2*max.subsets,.4*max.subsets,.6*max.subsets,.8*max.subsets,max.subsets),
         labels=c(0,20,40,60,80,100))
    if(sqrt.)
        label <- "normalized sqrt gcv or rss"
    else
        label <- "normalized gcv or rss"
    if(!zero.or.one.var)
        text(x=nrow(x) + 1.8, y=max.subsets/2, label,
             col=rh.col,
             xpd=NA, # no clip to plot region
             srt=90) # rotate text
    # bottom axis: variable names
    # axis() ignores the cex parameter (a bug?), so set cex globally, on.exit will restore it
    par(cex=cex.var)
    if(max.subsets == 0)
        axis(side=1, at=1, labels="intercept-only model")
    else
        axis(side=1, at=seq(1, nrow(x), by=1), labels=varlabs, las=3)
    invisible()
}
